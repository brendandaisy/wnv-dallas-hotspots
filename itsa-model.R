library(tidyverse)
library(lubridate)
library(sf)
library(terra)
library(tidyterra)
library(tidygeocoder)
library(cowplot)
library(usmap)

# ... expressions passed to filter
retry_coords <- function(df, method, ...) {
    retry_df <- filter(df, ...) |> 
        select(-long, -lat)
    
    df <- anti_join(df, retry_df)
    
    new_address_to_lat_longs <- geocode(
        retry_df, method=method, 
        street=address, city=city, postalcode=zip, country=country
    )
    
    combined_df <- bind_rows(df, new_address_to_lat_longs)
    
    return(combined_df)
}

pools0 <- readxl::read_xlsx("Comprehensive_WNV_DYS_tables_2024 (1).xlsx", sheet="Pools_updated") |> 
    transmute(
        date=as.Date(`Collection Date`),
        epiweek=epiweek(date),
        zone=Zone,
        address=Address, 
        city=City,
        zip=Zip,
        country="USA",
        num_tested=`# Tested`,
        result=case_match(
            Results,
            "N" ~ FALSE,
            "POS" ~ TRUE,
            .default=NA
        )
    )

distinct(pools0, Results)

pools0 |> 
    filter(is.na(result))

pools_geo1 <- geocode(pools0, street=address, city=city, postalcode=zip, country=country)

pools_geo2 <- retry_coords(
    pools_geo1, "arcgis", 
    long > -85 | long < -100 | is.na(lat) | is.na(long)
)

# TODO: there's definitely a little funkiness left still, with address, zip, and 
# geometry not always matching as expected. Probably related to a few duplicate addresses
# with different names
mpools_sf <- pools_geo2 |> 
    group_by(long, lat, month=month(date)) |>  # check that grouping by long lat ok here
    summarize(num_pos=sum(result), .groups="drop") |> 
    st_as_sf(coords=c("long", "lat"), crs="WGS84", remove=FALSE)

loc <- mpools_sf |> 
    distinct(long, lat) |> 
    st_as_sf(coords=c("long", "lat"), crs="WGS84")

mpools_sf |> 
    filter(month %in% 5:8) |> 
    ggplot() +
    geom_sf(aes(col=num_pos, size=num_pos), alpha=0.5) +
    facet_wrap(~month(month, label=TRUE, abbr=FALSE)) +
    scale_color_viridis_c() +
    scale_size(guide="none") +
    labs(col="num. positive\ntraps") +
    theme_bw()

ggsave("figs/summer-num-pos.pdf", width=6, height=5)

# pools_sf |> 
#     group_by(address, long, lat) |> 
#     summarise(any_pos=any(result), num_pos=sum(result)) |>
#     st_cast("POINT") |> 
#     ggplot() +
#     geom_sf(aes(col=any_pos, size=num_pos), alpha=0.5) +
#     theme_map()

# pools_geo2 |> 
#     distinct(long, lat) |> 
#     st_as_sf(coords=c("long", "lat"), crs="EPSG:3857") |> 
#     st_write("pools-coords.shp")

bioclim_list <- map(list.files("geofiles", "bio\\d+"), \(f) {
    rast(paste0("geofiles/", f)) |> 
        crop(st_bbox(loc), snap="out")
})
names(bioclim_list) <- paste0("bio", 1:19)

walk2(bioclim_list, 1:19, ~plot(.x, plg=list(title=paste0("bio", .y))))

# remove bio vars that had little to know variation within Dallas county:
bioclim_list <- bioclim_list[-c(10, 9, 2, 1)]

bioclim <- imap_dfc(bioclim_list, ~{
    tibble("{.y}":=extract(.x, mpools_sf)[,2])
}) |> 
    bind_cols(select(mpools_sf, long, lat)) |> 
    st_as_sf(coords=c("long", "lat"), crs="WGS84")

# quick check they seem to match
ggplot() +
    geom_spatraster(data=bioclim_list[[3]]) +
    geom_sf(aes(col=bio5), data=bioclim) +
    scale_color_viridis_c()

# haphazardly remove bios that are very correlated at the observed points
# TODO: may have to be redone if we get many more years of data
library(corrplot)
bc_corr <- cor(st_drop_geometry(bioclim))
corrplot(bc_corr)

# all the precipitation ones sound very similar, and indeed are correlated. bio15, 18, 19 seem
# the most interesting so keeping those
bioclim <- select(bioclim, -c(bio7, bio12, bio13, bio14, bio16, bio17))
bc_corr <- cor(st_drop_geometry(bioclim))
bc_corr < -0.8 | bc_corr > 0.8
corrplot(bc_corr)

bioclim_list <- discard_at(bioclim_list, c("bio7", "bio12", "bio13", "bio14", "bio16", "bio17"))

###
land_cover_labels <- function(lc) {
    lv <- c(11, 12, 21:24, 31, 41:43, 51, 52, 71:74, 81, 82, 90, 95) |> as.character()
    lab <- land_cover_names()
    ret <- factor(as.character(lc), levels=lv, labels=lab)
    return(fct_drop(ret))
}

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

lc_grid <- rast("geofiles/Annual_NLCD_LndCov_2024_CU_C1V1_mffklurgpwkbfm.tiff")
values(lc_grid) <- land_cover_labels(values(lc_grid))

lc <- mpools_sf |> 
    transmute(land_cover=extract(lc_grid, mpools_sf)[,2])

ggplot(lc) +
    geom_sf(aes(col=land_cover)) +
    theme_map()

# Add everything together!--------------------------------------------------------
mpools_obs <- mpools_sf |> 
    bind_cols(st_drop_geometry(bioclim), st_drop_geometry(lc))

# bioclim_grid <- project(bioclim_grid, "EPSG:3857")
covar_grid <- bioclim_list |> 
    rast() |> 
    as.data.frame(xy=TRUE) |> 
    as_tibble() |> 
    st_as_sf(coords=c("x", "y"), crs="WGS84")

covar_grid$land_cover <- extract(lc_grid, covar_grid)[,2]

ggplot(covar_grid) +
    geom_sf(aes(col=land_cover))

# combining the fit and prediction datasets:
# for non-spatial model:
# WARNING: covariates will now be centered and scaled. By definition this is 
# is relative to the included observation locations and times, and the prediction grid
# TODO: alternative is to define a "centering transformation" based on
# the covariate values at the observed traps/times and then apply that everywhere
mpools_join <- expand_grid(month=1:12, covar_grid) |> 
    mutate(num_pos=NA) |> 
    bind_rows(mpools_obs) |> 
    mutate(across(contains("bio"), ~(.x - mean(.x))/sd(.x)))

library(INLA)

mod <- as.formula(str_c(
    "num_pos ~ ",
    str_c(names(bioclim_list), collapse="+"),
    "+land_cover",
    "+f(month, model='rw2', cyclic=TRUE)"
))

fit <- inla(
    mod, data=st_drop_geometry(mpools_join),
    family="poisson",
    control.compute=list(dic=TRUE, return.marginals.predictor=TRUE),
    control.predictor=list(link=1),
    control.fixed=list(expand.factor.strategy="inla", prec=0.5, prec.intercept=0.01)
)
summary(fit)
fit$summary.random$month

pred_grid <- mpools_join |> 
    filter(is.na(num_pos)) |> 
    bind_cols(fit$summary.fitted.values[which(is.na(mpools_join$num_pos)),]) |> 
    st_as_sf()

# plot one "layer" (here, month) at a time:
pred_rast <- pred_grid |> 
    mutate(
        long=st_coordinates(pred_grid)[,1], 
        lat=st_coordinates(pred_grid)[,2]
    ) |>
    st_drop_geometry() |> 
    filter(month == 7) |> 
    select(long, lat, mean) |> 
    # pivot_wider(names_from=month, values_from=mean, names_prefix="mean_") |> 
    rast(crs="WGS84")

ggplot() +
    geom_spatraster(data=pred_rast) +
    geom_sf(aes(shape=num_pos > 0), filter(mpools_obs, month == 7), col="white") +
    scale_fill_viridis_c(option="inferno") +
    scale_shape_manual(values=c(1, 4)) +
    labs(shape="WNV detected", fill="predicted\npositive (mean)") +
    theme_map() +
    theme(legend.key=element_rect(fill="black"))

ggsave("figs/first-pred-map-july.png", width=5.4, height=4, bg="white")
