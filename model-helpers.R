# Build a prediction formula for predict()/generate() from a fitted bru object,
# summing every model component (fixed and random) on the linear predictor scale
# and wrapping it in the inverse link. Use `exclude` to drop components, e.g.
# exclude = "st_field" for covariate-only predictions.
build_pred_formula <- function(fit, exclude = NULL, link = "invlogit") {
    comps <- names(fit$bru_info$model$effects)
    bad   <- setdiff(exclude, comps)
    if (length(bad) > 0)
        stop("Not model components: ", paste(bad, collapse = ", "))
    
    eta <- paste(setdiff(comps, exclude), collapse = " + ")
    as.formula(paste0("~inla.link.", link, "(", eta, ")"), env = parent.frame())
}

# Build the non-separable spacetime field used as `st_field` in the model:
# a 1D temporal mesh over week_ind, a 2D spatial mesh over the bounding box of
# `fit_df`, and an INLAspacetime model "121" (critical diffusion) on them.
# Lengths are in km (project_crs()), times in weeks. Returns a list with the
# model object plus the meshes and boundary, which prediction grids need.
#
# TODO: a knot at every week made the temporal range too large; a better
# knot-spacing strategy is still needed. The coarse spatial mesh reduces
# compute time but biases the field; decrease edge lengths for production runs.
define_st_field <- function(
    fit_df,
    knot_spacing = 10,       # weeks between temporal knots
    edge_len     = 7,        # km, hexagon lattice spacing
    max_edge     = c(10, 10), # km
    prs          = c(10, 0.05), # P(spatial range < 10 km) = 0.05
    prt          = c(1,  0.05), # P(temporal range < 1 wk) = 0.05
    psigma       = c(1,  0.05)  # P(field sd > 1) = 0.05
) {
    wk_mesh <- fm_mesh_1d(seq(1, max(fit_df$week_ind), knot_spacing), degree = 1)

    bnd     <- st_as_sfc(st_bbox(fit_df), crs = st_crs(fit_df))
    loc_hex <- fm_hexagon_lattice(bnd, edge_len = edge_len)

    sp_mesh <- fm_mesh_2d(
        loc_hex, boundary = bnd,
        max.edge = max_edge,
        crs = fm_crs(fit_df)
    )

    st_model <- stModel.define(
        sp_mesh, wk_mesh, "121",
        control.priors = list(prs = prs, prt = prt, psigma = psigma),
        constr = TRUE
    )

    list(model = st_model, sp_mesh = sp_mesh, wk_mesh = wk_mesh, bnd = bnd)
}

# Find the `n_weeks` consecutive MMWR epiweeks (Sun-Sat, as in lubridate::epiweek,
# which matches the trap data) that best cover a calendar month: the window
# with the most days inside the month, ties broken by whichever window is
# centred closest to the middle of the month. Returns one row per day; weekends
# are dropped unless `include_weekends = TRUE`, since traps run Mon-Fri. The
# window is chosen on full weeks either way.
month_pred_dates <- function(month, year, n_weeks = 4, include_weekends = FALSE) {
    month_start <- make_date(year, month, 1)
    month_days  <- seq(month_start, ceiling_date(month_start, "month") - 1, by = "day")
    month_mid   <- mean(month_days)

    first_sun  <- floor_date(month_start, "week", week_start = 7)
    n_cand     <- ceiling((length(month_days) + 6) / 7) - n_weeks + 1
    cand_start <- first_sun + 7 * (seq_len(max(n_cand, 1)) - 1)

    window_days <- \(s) seq(s, by = "day", length.out = 7 * n_weeks)
    overlap     <- map_int(cand_start, \(s) sum(window_days(s) %in% month_days))
    mid_dist    <- map_dbl(cand_start, \(s) abs(as.numeric(mean(window_days(s)) - month_mid)))

    best <- cand_start[order(-overlap, mid_dist)[1]]

    dates <- window_days(best)
    if (!include_weekends)
        dates <- dates[wday(dates, week_start = 1) <= 5]

    tibble(date = dates) |>
        mutate(year = year(date), epiweek = epiweek(date))
}

# Cross every prediction pixel with every prediction date, and attach the
# model's temporal index `week_ind` by (year, epiweek) from `week_lookup`.
# Weeks absent from `week_lookup` get NA week_ind, with a warning, since the
# st_field has no time knot defined for them.
build_pred_df <- function(pix, pred_dates, week_lookup) {
    pred_dates <- left_join(pred_dates, week_lookup, by = c("year", "epiweek"))

    missing_wks <- unique(pred_dates$epiweek[is.na(pred_dates$week_ind)])
    if (length(missing_wks) > 0)
        warning("No week_ind in the fitted data for epiweek(s): ",
                paste(missing_wks, collapse = ", "))

    pix <- mutate(pix, pixel_id = row_number(), .before = 1)

    expand_grid(pred_dates, pix) |>
        relocate(pixel_id, week_ind, date, year, epiweek) |>
        st_as_sf()
}

# Average daily posterior draws within each pixel, separately for every draw.
# `draws` is the generate() matrix (one row per pixel-day, one column per
# draw); returns one row per pixel (rownames = pixel_id, sorted). The result
# is the posterior of the probability a trap placed on a random day of the
# week is positive, which keeps the full uncertainty for later summaries.
weekly_avg_draws <- function(draws, pixel_id) {
    n_days <- rowsum(rep(1, length(pixel_id)), pixel_id)
    rowsum(draws, pixel_id) / as.vector(n_days)
}

# Posterior summaries for each row of a draws matrix.
summarise_draws <- function(draws, probs = c(0.025, 0.5, 0.975)) {
    qs <- t(apply(draws, 1, quantile, probs = probs, names = FALSE))
    colnames(qs) <- paste0("q", probs)

    tibble(mean = rowMeans(draws), sd = apply(draws, 1, sd)) |>
        bind_cols(as_tibble(qs))
}

# ── Validation helpers ────────────────────────────────────────────────────────

# Reliability table: split predictions into `n_bins` equal-count bins and
# compare the mean predicted probability with the observed positive rate in
# each bin. obs_lo/obs_hi are Jeffreys 95% intervals for the observed rate.
calibration_bins <- function(p, y, n_bins = 10) {
    tibble(p, y) |>
        mutate(bin = ntile(p, n_bins)) |>
        group_by(bin) |>
        summarise(
            n      = n(),
            k      = sum(y),
            p_mean = mean(p),
            obs    = k / n,
            obs_lo = if_else(k == 0, 0, qbeta(0.025, k + 0.5, n - k + 0.5)),
            obs_hi = if_else(k == n, 1, qbeta(0.975, k + 0.5, n - k + 0.5)),
            .groups = "drop"
        )
}

# Area under the ROC curve: the probability a random positive (y = 1) gets a
# higher `p` than a random negative, computed from ranks (ties count half).
# NA when `y` has no positives or no negatives, where AUC is undefined.
auc <- function(p, y) {
    n_pos <- sum(y == 1)
    n_neg <- sum(y == 0)
    if (n_pos == 0 || n_neg == 0) return(NA_real_)
    (mean(rank(p)[y == 1]) - (n_pos + 1) / 2) / n_neg
}

# AUC of each posterior draw: `draws` is a generate() matrix with one row per
# observation and one column per draw. Returns one AUC per draw.
auc_draws <- function(draws, y) {
    apply(draws, 2, auc, y = y)
}

# Scores for binary predictions `p` of outcomes `y` (0/1). `p_ref` is a
# reference forecast (e.g. training-set positivity) for the Brier skill score:
# 1 = perfect, 0 = no better than the reference, < 0 = worse.
binary_scores <- function(p, y, p_ref) {
    tibble(
        n           = length(y),
        obs_rate    = mean(y),
        pred_rate   = mean(p),
        brier       = mean((p - y)^2),
        brier_skill = 1 - brier / mean((p_ref - y)^2),
        log_score   = -mean(y * log(p) + (1 - y) * log(1 - p)),
        auc         = auc(p, y)
    )
}
