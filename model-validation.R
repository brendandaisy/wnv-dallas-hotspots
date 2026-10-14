# model-comparison.R --------------------------------------------------------------
# Out-of-sample validation of the spatiotemporal WNV model in
# stf-cov_daymet_bioclim.R. Held-out predictions come from K-fold cross-validation over
# trap locations: each fold holds out every trap-day at ~1/K of the locations,
# so the score reflects predicting at an untrapped site, as the maps do.
# --------------------------------------------------------------------------------
library(INLA)
library(INLAspacetime)
library(inlabru)
library(tidyverse)
library(sf)
library(terra)
library(fmesher)
library(cowplot)
library(scoringRules)
library(ggthemes)

source("geo-helpers.R")
source("model-helpers.R")

# ── Data and model ─────────────────────────────────────────────────────────────
# Drops the few trap-days on Dec 31 of leap years, which Daymet doesn't cover.

fit_df <- read_sf("fit-df.gpkg") |>
    rename(geometry = geom) |>
    st_transform(project_crs()) |>
    drop_na(tmax, tmin, dayl, prcp) |>
    mutate(row_id=row_number(), before=1) |> 
    mutate(
        loc_id=vctrs::vec_group_id(as.data.frame(st_coordinates(geometry))), 
        .before = address
    )

st_field_def <- define_st_field(fit_df)
lc_grid      <- load_land_cover(target = vect(fit_df))

model <- result ~ Intercept(1, prec.linear = 0.5) +
    min_temp_cold_mon(bio06, model = "linear", prec.linear = 1) +
    precip_wet_mon(bio13, model = "linear", prec.linear = 1) +
    precip(prcp,  model = "linear", prec.linear = 1) +
    max_temp(tmax, model = "linear", prec.linear = 1) +
    land_cover(
        lc_grid, model = "iid",
        hyper = list(theta = list(prior = "pc.prec", param = c(1, 0.05)))
    ) +
    seasonal(
        epiweek, model = "rw2",
        cyclic=TRUE,
        values = 1:53, # all epiweeks, so weeks absent from the data still get an effect
        hyper = list(theta = list(prior = "pc.prec", param = c(1, 0.05)))
    ) +
    st_field(list(space = geometry, time = week_ind), model = st_field_def$model)

# - Top-level in-fit performance -
fit_all <- bru(
    model, fit_df, family = "binomial",
    options=bru_options(control.compute=list(dic=TRUE, mlik=FALSE, waic=TRUE, cpo=TRUE))
)

# make sure not too many failures for LPML
mean(fit_all$cpo$failure > 0)

tibble_row(
    `DIC (lower better)`=fit_all$dic$dic,
    `WAIC (lower better)`=fit_all$waic$waic,
    `LPML (higher better)`=sum(log(fit_all$cpo$cpo[fit_all$cpo$failure == 0]))
)

# CPO (leave-one-out prediction) patterns across space and time
fit_df_cpo <- fit_df |> 
    mutate(
        cpo=fit_all$cpo$cpo, 
        cpo_failure=fit_all$cpo$failure
    ) |> 
    filter(cpo_failure == 0)

ggplot(fit_df_cpo, aes(date, cpo)) +
    geom_point(alpha=0.1) +
    scale_x_date(date_labels="%b %Y", date_breaks="6 months") +
    labs(y="cpo (higher better)") +
    theme_bw()

fit_cpo_summ <- fit_df_cpo |> 
    group_by(loc_id, month=month(date, label=TRUE, abbr=TRUE)) |> 
    summarize(cpo=mean(cpo), .groups="drop")

ggplot(fit_cpo_summ) +
    geom_sf(aes(col=cpo), alpha=0.3) +
    facet_wrap(~month) +
    scale_colour_gradient2(low="tomato", mid="blue1", high="gray90", midpoint=0.5) +
    theme_bw()

# ── Cross-validation folds by trap location ───────────────────────────────────
# A location is a distinct trap coordinate; locations are randomly split into
# folds, so nearby traps can still land in the training set.
# All visits across time to an out-of-fold location are removed, matching a
# "fill in prediction map" task
# TODO discuss with group about this hold-out strategy for 1) model comparison
# and 2) quantifying forecast skill. We may want to try other things

n_folds <- 5

set.seed(2024)

fold_of_loc <- sample(rep_len(seq_len(n_folds), max(fit_df$loc_id)))
fit_df      <- mutate(fit_df, fold = fold_of_loc[loc_id])

count(st_drop_geometry(fit_df), fold, name = "trap_days")

# ── Held-out predictions ──────────────────────────────────────────────────────
# Refit on K-1 folds and predict the held-out trap-days. The posterior mean of
# p is the predictive probability that the trap is positive, which is what
# calibration assesses. Results are cached since each fold is a full refit.

cv_file <- "cv-trap-folds.rds"

if (file.exists(cv_file)) {
    cv_trap <- readRDS(cv_file)
} else {
    cv_trap <- map(sort(unique(fit_df$fold)), \(k) {
        train <- filter(fit_df, fold != k)
        test  <- filter(fit_df, fold == k)
        p_ref <- mean(train$result)
        
        fit_y <- bru(model, train, family = "binomial")
        
        draws <- generate(
            fit_y, test, build_pred_formula(fit_y),
            n.samples = 500, seed = 123+k
        )
        colnames(draws) <- str_c("draw_", 1:ncol(draws))
        
        mean_p_scores <- binary_scores(rowMeans(draws), test$result, p_ref) |>
            rename(mean_p_auc=auc)
        
        list(
            fold=k,
            mean_p_scores=mean_p_scores,
            auc_draws=tibble(
                draw = seq_len(ncol(draws)),
                auc  = auc_draws(draws, test$result)
            ),
            preds=draws |>
                as_tibble() |>
                mutate(row_id = test$row_id, y = test$result, .before = 1),
            mode_status=fit_y$mode$mode.status
        )
    })
    
    saveRDS(cv_trap, cv_file)
}

# Any fold whose fit didn't converge cleanly (should all be 0):
map_int(cv_trap, "mode_status")

# Shared layout for one score: held-out fold on y, optional no-skill line
score_panel <- function(
        data, y, title, subtitle = NULL, ref_line = NULL, ytitle=NULL) {
    p <- ggplot(data, aes(y = .data[[y]]))
    if (!is.null(ref_line))
        p <- p + geom_vline(xintercept = ref_line, linetype = "dashed",
                            col = "gray50", linewidth = 1.1)
    p +
        labs(title = title, subtitle = subtitle, x = NULL, y=ytitle) +
        theme_bw()
}

# Scores of the mean p: one value per year
mean_p_scores_trap <- map(cv_trap, "mean_p_scores") |>
    list_rbind(names_to = "fold")
mean_p_scores_trap

# Posterior median and intervals of per-draw AUC
auc_trap_summ <- map(cv_trap, "auc_draws") |>
    list_rbind(names_to="fold") |>
    group_by(fold) |>
    summarise(
        median   = median(auc),
        `q0.025` = quantile(auc, 0.025),
        `q0.25`  = quantile(auc, 0.25),
        `q0.75`  = quantile(auc, 0.75),
        `q0.975` = quantile(auc, 0.975)
    )

auc_trap_summ

# ── k-fold score panels ──────────────────────────────────────────────────────────

p_brier <- score_panel(
    mean_p_scores_trap, "fold",
    "Brier score", "mean p; lower is better",
    ytitle="fold"
) +
    geom_point(aes(x = brier), size = 1.6)

p_brier_skill <- score_panel(
    mean_p_scores_trap, "fold",
    "Brier skill score", 
    "mean p; baseline reference", 
    ref_line = 0, ytitle="fold"
) +
    geom_point(aes(x = brier_skill), size = 1.6)

p_log <- score_panel(
    mean_p_scores_trap, "fold",
    "log score", "mean p; lower is better",
    ytitle="fold"
) +
    geom_point(aes(x = log_score), size = 1.6)

p_auc <- score_panel(
    auc_trap_summ, "fold", 
    "AUC", "median, 50% and 95% over draws", 
    ref_line = 0.5, ytitle="fold"
) +
    geom_linerange(aes(xmin = `q0.025`, xmax = `q0.975`)) +
    geom_linerange(aes(xmin = `q0.25`,  xmax = `q0.75`), linewidth = 1.2) +
    geom_point(aes(x = median), size = 1.6)

plot_grid(
    p_brier, p_brier_skill, p_log, p_auc,
    ncol = 3, align = "hv"
)

# ══ One-week-ahead forecasting, 2024 ════════════════════════════════════════
# For every 2024 week with at least `min_week_traps` trap-days, refit on all
# data before that week (week_ind < target) and forecast every trap-day of the
# target week. Scores follow the leave-one-year-out section:
#
# - Brier, Brier skill and log score of the posterior mean p.
# - AUC of each of the 500 draws (NA in weeks with no positives).
# - CRPS of the week's positive count, from simulated trap-day outcomes.
#
# Reference forecast: the previous week's observed positivity (week_ind - 1,
# the most recent week with any trap data, even if it had fewer than
# `min_week_traps` trap-days and so isn't itself forecast), applied to every
# trap-day of the target week, and Binomial(trap-days, that positivity) for
# the weekly count.
#
# Per-week skill scores are unstable or undefined when the reference is near
# perfect (e.g. 0 positives in both weeks), so overall skill is pooled across
# forecast weeks: 1 - sum(model score) / sum(reference score).
#
# Note: the forecast week's covariates (tmax, prcp) are the observed values,
# so this assumes perfect knowledge of that week's weather.

fc_file        <- "1wk-forecast-2024.rds"
min_week_traps <- 20

fc_weeks <- st_drop_geometry(fit_df) |>
    filter(year == 2024) |>
    count(week_ind, epiweek) |>
    filter(n >= min_week_traps)

# The reference log score uses p_ref clamped away from 0 and 1, since the
# previous week's positivity can be exactly 0
clamp <- \(p) pmin(pmax(p, 1e-4), 1 - 1e-4)

if (file.exists(fc_file)) {
    cv_fc <- readRDS(fc_file)
} else {
    cv_fc <- map(fc_weeks$week_ind, \(wk) {
        train <- filter(fit_df, week_ind <  wk)
        test  <- filter(fit_df, week_ind == wk)
        p_ref <- mean(filter(fit_df, week_ind == wk - 1)$result)

        fit_w <- bru(model, train, family = "binomial")

        draws <- generate(
            fit_w, test, build_pred_formula(fit_w),
            n.samples = 500, seed = wk
        )
        colnames(draws) <- str_c("draw_", 1:ncol(draws))

        # Weekly positive count: simulate every trap-day outcome in each draw
        set.seed(wk)
        y_rep       <- matrix(rbinom(length(draws), 1, draws), nrow = nrow(draws))
        count_draws <- colSums(y_rep)
        ref_draws   <- rbinom(ncol(draws), nrow(test), p_ref)

        # All scores for the week, model and reference. Skill is -Inf when the
        # reference is perfect (0 positives in the previous and forecast week);
        # in those weeks the log skill depends on the clamp.
        y <- test$result

        week_scores <- binary_scores(rowMeans(draws), y, p_ref) |>
            rename(mean_p_auc = auc) |>
            mutate(
                week_ind      = wk,
                epiweek       = test$epiweek[1],
                obs           = sum(y),
                p_ref         = p_ref,
                .before       = 1
            ) |>
            mutate(
                brier_ref     = mean((p_ref - y)^2),
                log_score_ref = -mean(y * log(clamp(p_ref)) + (1 - y) * log(1 - clamp(p_ref))),
                log_skill     = 1 - log_score / log_score_ref,
                crps          = crps_sample(obs, matrix(count_draws, nrow = 1)),
                crps_ref      = crps_sample(obs, matrix(ref_draws,   nrow = 1)),
                crps_skill    = 1 - crps / crps_ref
            )

        list(
            week_ind    = wk,
            week_scores = week_scores,
            auc_draws   = tibble(
                draw = seq_len(ncol(draws)),
                auc  = auc_draws(draws, y)
            ),
            count_draws = count_draws,
            preds       = draws |>
                as_tibble() |>
                mutate(row_id = test$row_id, y = y, .before = 1),
            mode_status = fit_w$mode$mode.status
        )
    })

    saveRDS(cv_fc, fc_file)
}

names(cv_fc) <- map_chr(cv_fc, \(f) as.character(f$week_ind))

# Any forecast whose fit didn't converge cleanly (should all be 0):
table(map_int(cv_fc, "mode_status"))

# ── Forecast score summaries ──────────────────────────────────────────────────

# Per-week scores, model and reference
fc_week_scores <- map(cv_fc, "week_scores") |>
    list_rbind()

# Overall forecast skill, pooled over forecast weeks. Brier and log scores are
# per-trap-day means, so weekly values are weighted by trap-days (n), which
# equals pooling over all forecast trap-days. Each skill is computed before
# its score columns are overwritten by the pooled values.
fc_pooled <- fc_week_scores |>
    summarise(
        brier_skill   = 1 - sum(n * brier)     / sum(n * brier_ref),
        log_skill     = 1 - sum(n * log_score) / sum(n * log_score_ref),
        crps_skill    = 1 - sum(crps) / sum(crps_ref)
        # brier         = weighted.mean(brier,         n),
        # brier_ref     = weighted.mean(brier_ref,     n),
        # log_score     = weighted.mean(log_score,     n),
        # log_score_ref = weighted.mean(log_score_ref, n),
        # crps          = mean(crps),
        # crps_ref      = mean(crps_ref)
    )
fc_pooled

# Posterior median and intervals of per-draw AUC, for weeks with positives
auc_fc_summ <- map(cv_fc, "auc_draws") |>
    list_rbind(names_to = "week_ind") |>
    filter(!is.na(auc)) |>
    group_by(week_ind = as.integer(week_ind)) |>
    summarise(
        median   = median(auc),
        `q0.025` = quantile(auc, 0.025),
        `q0.25`  = quantile(auc, 0.25),
        `q0.75`  = quantile(auc, 0.75),
        `q0.975` = quantile(auc, 0.975)
    ) |>
    left_join(select(fc_week_scores, week_ind, epiweek), by = "week_ind")
auc_fc_summ

# Weekly count forecasts: predictive intervals vs observed and reference
fc_counts <- map(cv_fc, \(f) {
    qs <- quantile(f$count_draws, c(0.025, 0.25, 0.5, 0.75, 0.975))
    bind_cols(
        select(f$week_scores, week_ind, epiweek, n, obs, p_ref),
        tibble(`q0.025` = qs[1], `q0.25` = qs[2], median = qs[3],
               `q0.75` = qs[4], `q0.975` = qs[5])
    )
}) |>
    list_rbind()

# ── Forecast score panels ─────────────────────────────────────────────────────
# Weekly counts on top; per-week skill against the previous-week reference
# below, with the pooled skill in blue

# Shared layout for one score by forecast week; all panels span the same weeks
week_panel <- function(data, title, subtitle = NULL) {
    ggplot(data, aes(epiweek)) +
        scale_x_continuous(limits = range(fc_weeks$epiweek)) +
        labs(title = title, subtitle = subtitle, x = "2024 epiweek", y = NULL) +
        theme_bw()
}

fc_ref_col    <- "tomato"
fc_pooled_col <- "steelblue"
skill_floor   <- -2 # lower y limit of skill panels

# Per-week skill: values below `skill_floor` are drawn at the floor as
# triangles (open where -Inf, i.e. the reference was perfect); 0/0 weeks drop
skill_panel <- function(skill, title, pooled) {
    d <- fc_week_scores |>
        transmute(
            epiweek,
            skill = .data[[skill]],
            shown = pmax(skill, skill_floor),
            type  = case_when(
                skill == -Inf        ~ "reference perfect",
                skill <  skill_floor ~ "below axis",
                .default             = "skill"
            )
        ) |>
        filter(!is.nan(skill))

    week_panel(
        d, title,
        sprintf("pooled = %.2f; filled = below axis, open = reference perfect", pooled)
    ) +
        geom_hline(yintercept = 0, linetype = "dashed", col = "gray50", linewidth = 1.1) +
        geom_hline(yintercept = pooled, col = fc_pooled_col, linewidth = 0.8) +
        geom_line(aes(y = if_else(type == "skill", shown, NA)), col = "gray40",
                  na.rm = TRUE) +
        geom_point(aes(y = shown, shape = type), size = 1.6, fill = "black") +
        scale_shape_manual(
            values = c(skill = 16, "below axis" = 25, "reference perfect" = 6),
            guide = "none"
        ) +
        coord_cartesian(ylim = c(skill_floor, 1))
}

p_fc_counts <- week_panel(
    fc_counts, "Weekly positive count",
    "median, 50%, 95%; points = observed; red = reference"
) +
    geom_ribbon(aes(ymin = `q0.025`, ymax = `q0.975`), fill = "gray80") +
    geom_ribbon(aes(ymin = `q0.25`,  ymax = `q0.75`),  fill = "gray60") +
    geom_line(aes(y = median)) +
    geom_line(aes(y = n * p_ref), col = fc_ref_col, linetype = "dashed") +
    geom_point(aes(y = obs), size = 1.2)

p_fc_brier_skill <- skill_panel("brier_skill", "Brier skill score", fc_pooled$brier_skill)
p_fc_log_skill   <- skill_panel("log_skill",   "log score skill",   fc_pooled$log_skill)
p_fc_crps_skill  <- skill_panel("crps_skill",  "CRPS skill score",  fc_pooled$crps_skill)

p_fc_auc <- week_panel(
    auc_fc_summ, "AUC",
    "median, 50% and 95% over draws; weeks with positives"
) +
    geom_hline(yintercept = 0.5, linetype = "dashed", col = "gray50", linewidth = 1.1) +
    geom_linerange(aes(ymin = `q0.025`, ymax = `q0.975`)) +
    geom_linerange(aes(ymin = `q0.25`,  ymax = `q0.75`), linewidth = 1.2) +
    geom_point(aes(y = median), size = 1.6)

plot_grid(
    p_fc_counts,
    plot_grid(
        p_fc_brier_skill, p_fc_log_skill,
        p_fc_auc, p_fc_crps_skill,
        ncol = 2, align = "hv"
    ),
    ncol = 1, rel_heights = c(1, 2)
)

# ══ Leave-one-year-out: Bayesian scores ══════════════════════════════════════
# Zero-shot prediction of a whole year: each fold refits without one year and
# draws posterior samples of p for every trap-day of that year. Scores:
#
# - Brier, Brier skill and log score of the posterior mean p, which is the
#   predictive probability a trap is positive. Proper scores of single binary
#   outcomes only depend on this mean, so scoring each draw would just add a
#   penalty for posterior variance.
# - AUC of each draw, giving a posterior distribution of AUC. AUC only depends
#   on ranking within the year, so it ignores the overall positivity level.
# - CRPS of weekly positive counts during the WNV season (trap-days in
#   `wnv_months`). Each draw's trap-day outcomes are simulated and summed
#   within epiweek, giving a predictive distribution of the weekly count. Shared uncertainty (yearly level, seasonality, st_field) moves all
#   trap-days in a draw together, so this does reward well-calibrated
#   uncertainty. The CRPS skill reference is Binomial(trap-days in the week,
#   training positivity).
#
# For the held-out year the st_field is informed only through its temporal
# correlation with the neighbouring years, so this mostly tests the covariates
# and seasonal effect.

cv_file    <- "cv-year-folds.rds"
wnv_months <- 4:10 # Apr-Oct; weeks scored by CRPS

if (file.exists(cv_file)) {
    cv_year <- readRDS(cv_file)
} else {
    cv_year <- map(sort(unique(fit_df$year)), \(yr) {
        train <- filter(fit_df, year != yr)
        test  <- filter(fit_df, year == yr)
        p_ref <- mean(train$result)

        fit_y <- bru(model, train, family = "binomial")

        draws <- generate(
            fit_y, test, build_pred_formula(fit_y),
            n.samples = 500, seed = yr
        )
        colnames(draws) <- str_c("draw_", 1:ncol(draws))

        # Weekly positive counts in the WNV season: simulate every trap-day
        # outcome in each draw, then sum within epiweek (rows of count_draws
        # are sorted epiweeks)
        in_season <- month(test$date) %in% wnv_months

        set.seed(yr)
        y_rep       <- matrix(rbinom(length(draws), 1, draws), nrow = nrow(draws))
        count_draws <- rowsum(y_rep[in_season, ], test$epiweek[in_season])

        weekly <- st_drop_geometry(test) |>
            filter(in_season) |>
            group_by(epiweek) |>
            summarise(n = n(), obs = sum(result))

        ref_draws <- matrix(
            rbinom(nrow(weekly) * ncol(draws), weekly$n, p_ref),
            nrow = nrow(weekly)
        )

        weekly <- mutate(
            weekly,
            crps     = crps_sample(obs, count_draws),
            crps_ref = crps_sample(obs, ref_draws)
        )

        mean_p_scores <- binary_scores(rowMeans(draws), test$result, p_ref) |>
            rename(mean_p_auc=auc)

        list(
            mean_p_scores=mean_p_scores,
            auc_draws=tibble(
                draw = seq_len(ncol(draws)),
                auc  = auc_draws(draws, test$result)
            ),
            weekly_crps=weekly,
            preds=draws |>
                as_tibble() |>
                mutate(row_id = test$row_id, y = test$result, .before = 1),
            mode_status=fit_y$mode$mode.status
        )
    })
    names(cv_year) <- sort(unique(fit_df$year))

    saveRDS(cv_year, cv_file)
}

# Any fold whose fit didn't converge cleanly (should all be 0):
map_int(cv_year, "mode_status")

# ── Score summaries ───────────────────────────────────────────────────────────

# Scores of the mean p: one value per year
mean_p_scores_year <- map(cv_year, "mean_p_scores") |>
    list_rbind(names_to = "year")
mean_p_scores_year

# Posterior median and intervals of per-draw AUC
auc_year_summ <- map(cv_year, "auc_draws") |>
    list_rbind(names_to = "year") |>
    group_by(year) |>
    summarise(
        median   = median(auc),
        `q0.025` = quantile(auc, 0.025),
        `q0.25`  = quantile(auc, 0.25),
        `q0.75`  = quantile(auc, 0.75),
        `q0.975` = quantile(auc, 0.975)
    )

auc_year_summ

# CRPS per WNV-season week. Per-week skill is unstable in weeks where the
# reference CRPS is small, so the yearly summary pools: 1 - sum(crps) / sum(ref)
weekly_crps_year <- map(cv_year, "weekly_crps") |>
    list_rbind(names_to = "year") |>
    mutate(crps_skill = 1 - crps / crps_ref)

crps_year_summ <- weekly_crps_year |> 
    group_by(year) |> 
    reframe(
        crps_skill = 1 - sum(crps) / sum(crps_ref), # before crps is overwritten
        crps       = mean(crps)
    )

crps_year_summ

# ── Yearly score panels ──────────────────────────────────────────────────────────

p_brier <- score_panel(
    mean_p_scores_year, "year",
    "Brier score", "mean p; lower is better",
    ytitle = "held-out year"
) +
    geom_point(aes(x = brier), size = 1.6)

p_brier_skill <- score_panel(
    mean_p_scores_year, "year",
    "Brier skill score", "mean p; baseline reference",
    ref_line = 0, ytitle = "held-out year"
) +
    geom_point(aes(x = brier_skill), size = 1.6)

p_log <- score_panel(
    mean_p_scores_year, "year",
    "log score", "mean p; lower is better",
    ytitle = "held-out year"
) +
    geom_point(aes(x = log_score), size = 1.6)

p_auc <- score_panel(
    auc_year_summ, "year",
    "AUC", "median, 50% and 95% over draws",
    ref_line = 0.5, ytitle = "held-out year"
) +
    geom_linerange(aes(xmin = `q0.025`, xmax = `q0.975`)) +
    geom_linerange(aes(xmin = `q0.25`,  xmax = `q0.75`), linewidth = 1.2) +
    geom_point(aes(x = median), size = 1.6)

p_crps <- score_panel(
    weekly_crps_year, "year",
    "CRPS, weekly counts", "weeks Apr-Oct; lower is better",
    ytitle = "held-out year"
) +
    geom_boxplot(aes(x = crps), fill = "gray85", outliers = FALSE) +
    geom_jitter(aes(x = crps), height = 0.1, width = 0, size = 0.6, alpha = 0.5)

p_crps_skill <- score_panel(
    weekly_crps_year, "year",
    "CRPS skill score", "weeks Apr-Oct; baseline reference; x = pooled",
    ref_line = 0, ytitle = "held-out year"
) +
    geom_boxplot(aes(x = crps_skill), fill = "gray85", outliers = FALSE) +
    geom_jitter(aes(x = crps_skill), height = 0.1, width = 0, size = 0.6, alpha = 0.5) +
    geom_point(aes(x = crps_skill), crps_year_summ, shape = 4, size = 2, col = "tomato")

plot_grid(
    p_brier, p_brier_skill, p_log,
    p_auc, p_crps, p_crps_skill,
    ncol = 3, align = "hv"
)
