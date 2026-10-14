library(tidyverse)
library(sf)
library(terra)

# ── Project CRS ───────────────────────────────────────────────────────────────

# NAD83 / Texas North Central (EPSG:32138) with length unit km, so coordinates
# are exactly EPSG:32138 / 1000. Mesh edge lengths and the spatial range prior
# are in km. Note PROJ always reads x_0/y_0 in metres, whatever +units says.
project_crs <- function() {
    st_crs(paste(
        "+proj=lcc +lat_0=31.6666666666667 +lon_0=-98.5",
        "+lat_1=33.9666666666667 +lat_2=32.1333333333333",
        "+x_0=600000 +y_0=2000000 +datum=NAD83 +units=km +no_defs"
    ))
}

# ── Utilities ─────────────────────────────────────────────────────────────────

# Internal: tag each row with an integer `.loc_id` identifying its location.
# Rows sharing a geometry (e.g. one pixel across many dates) share an id, so
# rasters are extracted once per location and joined back on `.loc_id`.
.add_loc_id <- function(points) {
    points <- st_as_sf(points)
    xy     <- as.data.frame(st_coordinates(points))
    mutate(points, .loc_id = vctrs::vec_group_id(xy))
}

# Internal: one row per distinct location, carrying only `.loc_id` and geometry.
.distinct_locs <- function(points) {
    points |>
        select(.loc_id) |>
        distinct(.loc_id, .keep_all = TRUE)
}

# Internal: points as a SpatVector in the CRS of raster `r`, so values are read
# from the raster's native grid instead of a resampled copy.
.points_in_crs <- function(points, r) {
    project(vect(points), crs(r))
}

# ── Land cover (NLCD) ─────────────────────────────────────────────────────────

land_cover_names <- function() {
    c(
        "Open Water",
        "Perennial Ice/Snow",
        "Developed, Open Space",
        "Developed, Low Intensity",
        "Developed, Medium Intensity",
        "Developed, High Intensity",
        "Barren Land",
        "Deciduous Forest",
        "Evergreen Forest",
        "Mixed Forest",
        "Dwarf Scrub",
        "Shrub/Scrub",
        "Grassland/Herbaceous",
        "Sedge/Herbaceous",
        "Lichens",
        "Moss",
        "Pasture/Hay",
        "Cultivated Crops",
        "Woody Wetlands",
        "Emergent Herbaceous Wetlands"
    )
}

land_cover_labels <- function(lc, drop = FALSE) {
    ret <- recode_values(
        as.character(lc),
        from = as.character(c(11, 12, 21:24, 31, 41:43, 51, 52, 71:74, 81, 82, 90, 95)),
        to   = land_cover_names()
    ) |>
        factor()

    if (drop) fct_drop(ret) else ret
}

#' Load the NLCD land cover raster, relabel values, and reproject to match target
#'
#' Reprojection uses nearest-neighbour so land cover classes are never averaged.
#'
#' @param path    Path to the NLCD .tiff file.
#' @param target  An sf, SpatVector, or SpatRaster whose CRS is used as target.
#' @return SpatRaster with labelled factor values reprojected to target CRS.
load_land_cover <- function(
    path = "geofiles/Annual_NLCD_LndCov_2024_CU_C1V1_mffklurgpwkbfm.tiff",
    target,
    drop=FALSE
) {
    lc <- rast(path)
    values(lc) <- land_cover_labels(values(lc), drop = drop)
    project(lc, crs(target), method = "near")
}


# ── Bioclim (annual PRISM-style rasters) ──────────────────────────────────────

# Internal: extract one bioclim raster at point locations, tagging rows with year.
.extract_bioclim_points <- function(path, year, points) {
    r    <- rast(path)
    vals <- terra::extract(r, .points_in_crs(points, r), ID = FALSE)
    st_as_sf(points) |>
        mutate(rvals = as.double(unlist(vals)), year = as.integer(year))
}

#' Load one bioclim variable for one or more years at point locations
#'
#' Returns a long-format sf suitable for joining onto other data by
#' (location, year). Points may be in any CRS; they are transformed to the
#' raster's CRS for extraction.
#'
#' File naming convention: `<var>_equiv_20<YY>.tiff` under geofiles_dir.
#'
#' @param var         Bioclim variable name, e.g. "bio05".
#' @param years       4-digit year(s) as integer or character, e.g. c(2021, 2022).
#' @param points      sf with point locations. Will be passed as-is to
#'                    terra::extract(); use .distinct_locs() upstream for efficiency.
#' @param geofiles_dir Directory containing the *_equiv_*.tiff files.
#' @return sf data frame with all columns of `points` plus `year` and `<var>`.
load_bioclim_var <- function(var, years, points, geofiles_dir = "geofiles") {
    years   <- as.integer(years)
    yr2     <- formatC(years %% 100, width = 2, flag = "0")
    pattern <- sprintf("^%s_equiv_20(%s)\\.tiff$", var, paste(yr2, collapse = "|"))
    files   <- sort(list.files(geofiles_dir, pattern = pattern, full.names = TRUE))

    if (length(files) == 0)
        stop(sprintf("No bioclim files found for '%s' in '%s'", var, geofiles_dir))

    varname <- gsub("-", "_", var)

    map2_dfr(files, sort(years), .extract_bioclim_points, points = points) |>
        rename(!!varname := rvals)
}

#' Attach multiple bioclim variables to a point sf data frame
#'
#' Works for both trap observations (many rows per location, multi-year) and
#' prediction pixel grids (one row per location, typically single year).
#' Extractions are performed at distinct spatial locations for efficiency, then
#' joined back by (location, year) when `points` has a year column, or just
#' by location otherwise (appropriate when a single year is specified).
#'
#' @param points      sf data frame: trap observations or prediction pixels.
#' @param biovars     Character vector of bioclim variable names.
#' @param years       4-digit year(s). When `points` has a `year` column all
#'                    years are loaded and matched; otherwise a single year
#'                    should be given and values are joined by location alone.
#' @param filename    Optional file path for caching the result (e.g. ".gpkg",
#'                    ".shp"). Format is inferred from the extension.
#' @param save        If TRUE and filename is not NULL, write the result to disk.
#' @param geofiles_dir Directory containing the *_equiv_*.tiff files.
#' @return `points` with one new column per bioclim variable.
add_bioclim_covariates <- function(
    points,
    biovars,
    years,
    filename     = NULL,
    save         = TRUE,
    geofiles_dir = "geofiles"
) {

    points        <- .add_loc_id(points)
    distinct_locs <- .distinct_locs(points)
    has_year_col  <- "year" %in% names(points)

    for (bvar in biovars) {
        bio_df  <- load_bioclim_var(bvar, years, distinct_locs, geofiles_dir)
        varname <- gsub("-", "_", bvar)

        if (has_year_col) {
            # Multi-year trap data: join by location + year
            join_df  <- st_drop_geometry(bio_df)
            join_key <- c(".loc_id", "year")
        } else {
            # Single-year prediction pixels: drop year col, join by location only
            join_df  <- select(st_drop_geometry(bio_df), -year)
            join_key <- ".loc_id"
        }

        points <- left_join(points, join_df, by = join_key)
    }

    points <- select(points, -.loc_id)

    if (save && !is.null(filename))
        write_sf(points, filename)

    points
}


# ── Daymet (daily rasters) ────────────────────────────────────────────────────

#' Load a multi-layer Daymet raster for one variable and year
#'
#' @param var         Variable name matching file prefix, e.g. "tmax".
#' @param year_suffix 2-digit year suffix used in filenames, e.g. 21 for 2021.
#' @param daymet_dir  Directory containing the Daymet .tif files.
#' @return Multi-layer SpatRaster (one layer per day-of-year).
load_daymet_raster <- function(var, year_suffix,
                                daymet_dir = "geofiles/Daymet/Daymet_Daily") {
    rast(file.path(daymet_dir, sprintf("%s_%s.tif", var, year_suffix)))
}

# Internal: extract daily Daymet values at distinct locations, returning a
# long-format data frame with columns (.loc_id, date, <var>).
.extract_daymet_long <- function(r, points, var) {
    names(r) <- as.character(as.Date(time(r)))

    vals <- terra::extract(r, .points_in_crs(points, r), ID = FALSE)

    bind_cols(st_drop_geometry(points) |> select(.loc_id), vals) |>
        pivot_longer(
            cols      = -.loc_id,
            names_to  = "date",
            values_to = var
        ) |>
        mutate(date = as.Date(date))
}

#' Attach multiple Daymet variables to a point sf data frame
#'
#' Works for both trap observations and prediction pixel grids, provided
#' `points` has a `date` column. Extractions are performed at distinct spatial
#' locations for efficiency, then joined back by (location, date).
#'
#' @param points        sf data frame with a `date` column.
#' @param dayvars       Character vector of Daymet variable names,
#'                      e.g. c("tmax", "tmin", "dayl", "prcp").
#' @param years         Four digit years
#' @param filename      Optional file path for caching the result.
#' @param save          If TRUE and filename is not NULL, write the result to disk.
#' @param daymet_dir    Directory containing the Daymet .tif files.
#' @return `points` with one new column per Daymet variable.
add_daymet_covariates <- function(
    points,
    dayvars,
    years,
    filename   = NULL,
    save       = TRUE,
    daymet_dir = "geofiles/Daymet/Daymet_Daily"
) {

    year_suffixes <- years %% 100
    points        <- .add_loc_id(points)
    distinct_locs <- .distinct_locs(points)

    for (var in dayvars) {
        long_df <- map_dfr(year_suffixes, function(yr) {
            r <- load_daymet_raster(var, yr, daymet_dir)
            .extract_daymet_long(r, distinct_locs, var)
        })
        points <- left_join(points, long_df, by = c(".loc_id", "date"))
    }

    points <- select(points, -.loc_id)

    if (save && !is.null(filename))
        write_sf(points, filename)

    points
}
