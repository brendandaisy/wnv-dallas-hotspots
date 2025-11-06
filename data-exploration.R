library(tidyverse)
library(lubridate)
library(sf)
library(terra)
library(tidyterra)
library(tidygeocoder)
library(cowplot)

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
    ) |> 
    # fix the presumably erroneous 2022 dates:
    mutate(date=as.Date(str_replace(date, "2022", "2024")))

pools0 |> 
    filter(is.na(result))

pools_geo1 <- geocode(pools0, street=address, city=city, postalcode=zip, country=country)

pools_geo2 <- retry_coords(
    pools_geo1, "arcgis", 
    long > -85 | long < -100 | is.na(lat) | is.na(long)
)

pools_sf <- pools_geo2 |> 
    st_as_sf(coords=c("long", "lat"), crs="EPSG:5070", remove=FALSE) |> 
    st_transform("EPSG:5070")

loc <- pools_sf |> 
    distinct(long, lat, geometry)

# temporal trend - identifying a useful unit of analysis--------------------------
library(zoo)

pools_date <- pools_sf |> 
    st_drop_geometry() |> 
    group_by(date) |> 
    summarize(num_pos=sum(result)) |> 
    mutate(pos_wk=rollapply(num_pos, 7, mean, fill=NA))

ggplot(pools_date, aes(date)) +
    geom_point(aes(y=num_pos), alpha=0.2) +
    geom_line(aes(y=pos_wk), col="blue") +
    theme_half_open()

## number of visits each week
pools_sf |> 
    st_drop_geometry() |> 
    count(epiweek) |> 
    ggplot(aes(epiweek, n)) +
    geom_col() +
    theme_half_open()

## variability within a given week and site
wpools_sf <- pools_sf |> 
    group_by(long, lat, epiweek, year=year(date)) |> 
    summarize(num_visits=n(), num_pos=sum(result), .groups="drop") |> 
    mutate(any_pos=ifelse(num_pos > 0, 1, 0), prop=num_pos / num_visits)

# "just under 5% of sites were visited multiple times in a week"
nrow(filter(wpools_sf, num_visits > 1)) / nrow(wpools_sf)

# "of these, 17% had both positive and negative results in a week"
nrow(filter(wpools_sf, prop > 0, prop < 1)) / nrow(filter(wpools_sf, num_visits > 1))

# (based on this, choosing to go with weekly resolution and a binary response!)

#  spatiotemporal variogram-------------------------------------------------------
# (TODO?)

# total number of positives per month per site----------------------------------
mpools_sf <- pools_sf |> 
    group_by(long, lat, month=month(date)) |>  # check that grouping by long lat ok here
    summarize(num_visits=n(), num_pos=sum(result), .groups="drop")

mpools_sf |> 
    filter(month %in% 3:11) |> 
    ggplot() +
    geom_sf(aes(size=num_visits, col=ifelse(num_pos == 0, NA, num_pos)), alpha=0.6) +
    facet_wrap(~month(month, label=TRUE, abbr=FALSE)) +
    scale_color_viridis_c() +
    # scale_size(guide="none") +
    labs(size="num. visits", col="num. positive\ntraps") +
    theme_bw()

# there's a modest association between num. positive and num. visits, so we
# will likely want to change the response variable in some way or control for this
ggplot(mpools_sf, aes(num_visits, num_pos)) +
    geom_point(alpha=0.2) +
    geom_smooth(method="lm", se=FALSE, col="blue")

# Do same sites tend to stay positive between weeks during WNV season?------------
library(tseries)

site_stats <- wpools_sf |> 
    filter(epiweek > epiweek("2024-06-15"), epiweek < epiweek("2024-09-01")) |> 
    group_by(long, lat) |> 
    reframe(
        num_visits=sum(num_visits),
        prop=sum(any_pos)/num_visits,
        ac=cor(any_pos, lag(any_pos), use="pairwise.complete.obs"),
        max_run=max(rle(any_pos)$lengths[rle(any_pos)$values == 1])
    )

max_vis <- max(site_stats$num_visits)

site_stats |> 
    arrange(prop) |> 
    mutate(site=fct_inorder(interaction(long, lat))) |> 
    ggplot(aes(site, prop, group=1)) +
    geom_point(aes(y=num_visits/max_vis), col="blue", alpha=0.3) +
    geom_line() +
    geom_point() +
    scale_y_continuous(sec.axis=sec_axis(transform=~.*max_vis, name="num. visits")) +
    theme_half_open() +
    theme(axis.text.x=element_blank(), axis.ticks.x=element_blank())

site_stats |> 
    filter(!is.na(ac)) |> 
    ggplot(aes(ac)) +
    geom_histogram() +
    theme_half_open()

# spatiotemporal "base" model-----------------------------------------------------
library(INLA)
library(INLAspacetime)
library(fmesher)


plot(loc_grid)
st_distance(loc_grid)

## define the temporal mesh
wk_mesh <- fm_mesh_1d(1:52, degree=1)
wk_mesh$n

## define the spatial mesh
bnd <- st_as_sfc(st_bbox(loc), crs=st_crs(loc))
loc_hex <- fm_hexagon_lattice(bnd, edge_len=0.03)

sp_mesh <- fm_mesh_2d(
    loc_hex, boundary=bnd, 
    max.edge=c(0.03, 0.1), 
    crs=fm_crs(loc)
)

ggplot() +
    geom_fm(data=sp_mesh) +
    geom_sf(data=loc)

matern <- inla.spde2.pcmatern(sp_mesh, prior.sigma=c(3, 0.01), prior.range=c(0.5, 0.01))

mod <- any_pos ~ field(geometry, model=matern) + Intercept(1)

# what was a week with a lot of positives?
wpools_sf |> 
    group_by(epiweek, year) |> 
    reframe(num_pos=sum(num_pos)) |> 
    ggplot(aes(epiweek)) +
    geom_point(aes(y=num_pos)) +
    theme_half_open()

quick_fit_epiweek <- function(epiweek) {
    wpools_sub <- filter(wpools_sf, epiweek == {{epiweek}})
    fit <- bru(mod, wpools_sub, family="binomial")
    
    pix <- fm_pixels(sp_mesh, dims=c(150, 150), mask=bnd)
    pred <- predict(fit, pix, ~inla.link.invlogit(field + Intercept), n.samples=1000)
    list(ew=epiweek, data=wpools_sub, pred=pred)
}

plot_quick_fit <- function(ff) {
    ggplot() +
        geom_sf(data=bnd) +
        gg(ff$pred, aes(fill=mean), geom="tile") +
        geom_sf(aes(shape=as.factor(any_pos)), ff$data, col="gray70") +
        scale_fill_viridis_c(option="inferno") +
        scale_shape_manual(values=c(1, 4)) +
        labs(title=paste0("epiweek ", ff$ew), fill="prediction\nmean", shape="WNV positive") +
        theme_map()
}

f1 <- quick_fit_epiweek(28)
f2 <- quick_fit_epiweek(35)
f3 <- quick_fit_epiweek(40)

plots <- map(list(f1, f2, f3), plot_quick_fit)
plot_grid(plotlist=plots, nrow=1, align="h")

## older plotting stuff
pred <- predict(fit, pix, ~field + Intercept, n.samples=1000)

p1 <- ggplot() +
    geom_sf(data=bnd) +
    gg(pred, aes(fill=mean), geom="tile") +
    geom_sf(aes(shape=as.factor(any_pos)), wpools_sub, col="gray50") +
    scale_fill_viridis_c(option="inferno") +
    scale_shape_manual(values=c(1, 4)) +
    labs(title="epiweek 40", fill="prediction\nmean", shape="WNV positive") +
    theme_map()
    # theme(legend.key=element_rect(fill="gray30", color=NA))

p2 <- ggplot() +
    geom_sf(data=bnd) +
    gg(pred, aes(fill=sd), geom="tile") +
    geom_sf(aes(shape=as.factor(any_pos)), wpools_sub, col="gray50") +
    scale_fill_viridis_c(option="mako") +
    scale_shape_manual(values=c(1, 4)) +
    labs(fill="prediction\nstd.", shape="WNV positive") +
    theme_map() +
    theme(legend.key=element_rect(fill="gray30", color=NA))

plot_grid(p1, p2, nrow=1, align="h")
