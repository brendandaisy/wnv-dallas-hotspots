# stf-cov_daymet_bioclim.R -------------------------------------------------------
# Fit a non-separable spatiotemporal process for WNV vector index using
# INLAbru / INLAspacetime, with Daymet daily and Bioclim annual covariates.
# --------------------------------------------------------------------------------
library(INLA)
library(INLAspacetime)
library(inlabru)
library(tidyverse)
library(lubridate)
library(sf)
library(terra)
library(tidyterra)
library(tidygeocoder)
library(cowplot)
library(furrr)
library(nngeo)
library(fmesher)
library(ncdf4)

source("geo-helpers.R")
source("model-helpers.R")

# ── Load trap observations ─────────────────────────────────────────────────────
# Everything downstream (mesh, priors, prediction grid) works in the project
# CRS: Texas North Central in km. The lat/long columns stay as WGS 84 degrees.

traps_sf <- read_sf("traps-all-years.shp") |>
    st_transform(project_crs()) |>
    rename(num_trapped = nm_trpp) |>
    mutate(year = year(date), .after = epiweek) |>
    group_by(year, epiweek) |>
    mutate(week_ind = cur_group_id()) |>
    ungroup()

# ── Bioclim covariates (annual) ────────────────────────────────────────────────
# Cached to traps_sf_bioclim.gpkg; re-extraction is skipped on subsequent runs.
covyears <- c(2021, 2022, 2023, 2024)

biovars  <- c(
    "bio05", "bio06", "bio08", "bio10", "bio11", 
    "bio13", "bio14", "bio18", "bio19"
)

fit_df <- add_bioclim_covariates(
    traps_sf, biovars, covyears,
    save=FALSE
)

# ── Daymet covariates (daily) ──────────────────────────────────────────────────
dayvars <- c("tmax", "tmin", "dayl", "prcp")

fit_df <- add_daymet_covariates(
    fit_df, dayvars, covyears,
    save=TRUE,
    filename="fit-df.gpkg"
)

# instead of extracting covariates, can load preexisting fit_df directly
# (st_transform guards against a cache written before the switch to km):
# drop a few 2024-12-31 observations with NA daymet values from leap year
fit_df <- read_sf("fit-df.gpkg") |> 
    rename(geometry=geom) |>
    st_transform(project_crs()) |> 
    drop_na(tmax, tmin, dayl, prcp)

# ── Quick covariate diagnostic plots ──────────────────────────────────────────

ggplot(fit_df, aes(year, bio19, col = as.factor(result))) +
    geom_boxplot()

ggplot(fit_df, aes(epiweek, prcp, col = as.factor(result))) +
    geom_point(alpha = 0.5, size = 1.05)

ggplot(fit_df, aes(epiweek, tmax, col = as.factor(result))) +
    geom_point(alpha = 0.5, size = 1.05)

# ── Spacetime field ────────────────────────────────────────────────────────────
# Meshes and the INLAspacetime model "121" (see define_st_field()).
st_field_def <- define_st_field(fit_df)

sp_mesh <- st_field_def$sp_mesh
bnd     <- st_field_def$bnd

# Visually check: triangles should be roughly even, with a decent buffer
# around observation points, and small enough to capture local variation.
ggplot() +
    geom_fm(data = sp_mesh) +
    geom_sf(data = fit_df)

# ── Land cover (NLCD) ─────────────────────────────────────────────────────────

lc_grid <- load_land_cover(target = vect(fit_df), drop=TRUE)
plot(lc_grid)

# ── Model formula and fit ──────────────────────────────────────────────────────

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

fit <- bru(model, fit_df, family = "binomial")

# sometimes there are divergence warnings. These aren't great, but its
# definitely bad if this is something other than 0:
fit$mode$mode.status

summary(fit)

# plot the posterior seasonal effect, if present:
if (!is.null(fit$summary.random$seasonal))
    plot(fit, "seasonal")

# ── Covariate coefficient forest plot ─────────────────────────────────────────

lc_visited_labs <- levels(lc_grid)[[1]]$label
names(fit$marginals.random$land_cover) <- lc_visited_labs

covariates_marg <- c(fit$marginals.fixed, fit$marginals.random$land_cover)[-1]

covariates_summ <- map(covariates_marg, \(m) {
    hdpi <- inla.hpdmarginal(c(0.5, 0.95), m) |>
        as_tibble() |>
        mutate(p = c(0.5, 0.95)) |>
        pivot_wider(names_from = p, values_from = c(low, high))
    mutate(hdpi, mean = inla.emarginal(\(x) x, m))
}) |>
    list_rbind(names_to = "var")

ggplot(covariates_summ, aes(y = var)) +
    geom_vline(xintercept = 0, linetype = "dashed", col = "gray50", linewidth = 1.1) +
    geom_linerange(aes(xmin = `low_0.95`, xmax = `high_0.95`)) +
    geom_linerange(aes(xmin = `low_0.5`,  xmax = `high_0.5`), linewidth = 1.2) +
    geom_point(aes(x = mean), size = 1.6) +
    labs(x = "log odds", y = NULL)

# ── Weekly prediction maps ─────────────────────────────────────────────────────

pred_month <- 7
pred_year  <- 2024

pix_pred <- fm_pixels(sp_mesh, dims = c(100, 100), mask = bnd)

# Four MMWR weeks (Sun-Sat) best covering the month, one row per weekday
pred_dates <- month_pred_dates(pred_month, pred_year, n_weeks = 4, include_weekends=FALSE)

# Map (year, epiweek) to the week_ind used by st_field in the fit
week_lookup <- fit_df |>
    st_drop_geometry() |>
    distinct(year, epiweek, week_ind)

# One row per pixel x day
pred_df <- build_pred_df(pix_pred, pred_dates, week_lookup)

# Extract Bioclim and Daymet precipitation for prediction pixels -------------------

pred_df <- add_bioclim_covariates(
    pred_df, biovars, pred_year,
    save=FALSE
)

pred_df <- add_daymet_covariates(
    pred_df, dayvars, pred_year,
    save=FALSE
    # save=TRUE,
    # filename="pred-df.gpkg"
)

# should all be the same
ggplot() +
    gg(pred_df, aes(fill=bio06), geom="tile") +
    facet_wrap(~date)
   
# should be ~different each date 
ggplot() +
    gg(pred_df, aes(fill=prcp), geom="tile") +
    facet_wrap(~date)

# Generate posterior predictive samples ---------------------------------------
pred_formula <- build_pred_formula(fit)

# Run generate() one week at a time to limit memory, reducing each week's daily
# draws to weekly-average draws (pixels x samples) before the next week. The
# fixed seed gives the same posterior draws in every call, so sample j is one
# joint draw across weeks.
pred_seed <- 123
weekly_draws <- split(pred_df, pred_df$epiweek) |>
    map(\(wk_df) {
        draws <- generate(
            fit,
            wk_df,
            pred_formula,
            n.samples = 1000,
            seed = pred_seed
        )
        weekly_avg_draws(draws, wk_df$pixel_id)
    })

# Posterior summaries per pixel-week, joined back onto the pixel geometry
pix_pred_id <- mutate(pix_pred, pixel_id = row_number(), .before = 1)

pred_weekly <- weekly_draws |>
    imap(\(d, wk) summarise_draws(d) |>
        mutate(pixel_id = as.integer(rownames(d)), epiweek = as.integer(wk),
               .before = 1)) |>
    list_rbind() |>
    left_join(pix_pred_id, by = "pixel_id") |>
    st_as_sf()

# Posterior mean probability a trap set on a weekday is positive
# TODO update all maps to be greyed out if there are bodies of water etc
ggplot() +
    gg(pred_weekly, aes(fill = mean), geom = "tile") +
    geom_sf(
        aes(shape=ifelse(result == 1, "positive", "negative")),
        filter(traps_sf, year == pred_year, epiweek %in% unique(pred_df$epiweek)),
        col="gray70"
    ) +
    facet_wrap(~epiweek, labeller=label_both) +
    scale_fill_viridis_c(name = "predicted") +
    scale_shape_manual(values=c(1, 4), name="observed") +
    scale_x_continuous(expand=expansion()) +
    scale_y_continuous(expand=expansion()) +
    labs(
        subtitle=paste0(
            "Mean predicted prevalence for ",
            month(pred_month, label=TRUE, abbr=FALSE),
            ", ",
            pred_year
        )
    ) +
    theme_bw() +
    theme(
        axis.title=element_blank(),
        axis.text=element_blank(),
        axis.ticks=element_blank()
    )

# Exceeding probabilities above a specified threshold
# Threshold: observed positivity of all trap-days in pred_month of pred_year,
# so the map shows hotspots relative to that month.
# P(exceed) is the share of posterior draws whose weekly-average
# probability is above it.
pos_threshold <- traps_sf |>
    filter(month(date) == pred_month, year == pred_year) |>
    pull(result) |>
    mean()

exceed_weekly <- weekly_draws |>
    imap(\(d, wk) tibble(
        pixel_id = as.integer(rownames(d)),
        epiweek  = as.integer(wk),
        p_exceed = rowMeans(d > pos_threshold)
    )) |>
    list_rbind()

pred_weekly <- pred_weekly |>
    select(-any_of("p_exceed")) |>
    left_join(exceed_weekly, by = c("pixel_id", "epiweek"))

ggplot() +
    gg(pred_weekly, aes(fill = p_exceed), geom = "tile") +
    geom_sf(
        aes(shape=ifelse(result == 1, "positive", "negative")),
        filter(traps_sf, year == pred_year, epiweek %in% unique(pred_df$epiweek)),
        col="gray70"
    ) +
    facet_wrap(~epiweek, labeller=label_both) +
    scale_fill_viridis_c(
        name = "exceedance\nprobability",
        option = "magma", limits = c(0, 1)
    ) +
    scale_shape_manual(values=c(1, 4), name="observed") +
    labs(
        subtitle=paste0(
            "Probability of weekly risk exceeding ",
            round(100 * pos_threshold, 1),
            "% threshold for ",
            month(pred_month, label=TRUE, abbr=FALSE),
            ", ",
            pred_year
        )
    ) +
    scale_x_continuous(expand=expansion()) +
    scale_y_continuous(expand=expansion()) +
    theme_bw() +
    theme(
        axis.title=element_blank(),
        axis.text=element_blank(),
        axis.ticks=element_blank()
    )

# ── Spatiotemporal hyperparameter posteriors ───────────────────────────────────
post1 <- inla.tmarginal(
    \(x) exp(x), 
    fit$internal.marginals.hyperpar$`Theta1 for st_field`
) |>
    as_tibble() |> mutate(parameter = "spatial range (km)")

post2 <- inla.tmarginal(
    \(x) exp(x), 
    fit$internal.marginals.hyperpar$`Theta2 for st_field`
) |>
    as_tibble() |> mutate(parameter = "temporal range (wk)")

post3 <- inla.tmarginal(
    \(x) exp(x), 
    fit$internal.marginals.hyperpar$`Theta3 for st_field`
) |>
    as_tibble() |> mutate(parameter = "field sd.")

bind_rows(post1, post2, post3) |>
    mutate(parameter = fct_inorder(parameter)) |>
    ggplot(aes(x, y, col = parameter)) +
    geom_line() +
    labs(x=NULL, y="posterior density") +
    facet_wrap(~parameter, scales = "free") +
    theme_half_open() +
    theme(
        legend.position="none"
    )
