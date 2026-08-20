# spatiotemporal-fit.R------------------------------------------------------------
# fit a baseline non-separable spatiotemporal process using INLAspacetime---------
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
library(dplyr)
library(raster, exclude=c("select"))
library(furrr)
library(nngeo)
library(fmesher)
library(ncdf4)

setwd("~/WNV")

source("geo-helpers.R")

traps_sf <- read_sf("traps-all-years.shp") |> 
  rename(num_trapped=nm_trpp) |> 
  mutate(year=year(date), .after=epiweek) |> 
  group_by(year, epiweek) |> 
  mutate(week_ind=cur_group_id()) |> 
  ungroup()
length(unique(traps_sf$geometry))


loc <- distinct(traps_sf, lat, long, geometry)
class(loc)
wk <- seq(min(traps_sf$week_ind), max(traps_sf$week_ind), 1)

## define the temporal mesh
# TODO no idea what a "good" choice for knots is, but a knot at every week
# made the temporal range way to big. Overall, I don't really understand the
# temporal behavior of the model very well yet
wk_mesh <- fm_mesh_1d(seq(1, max(wk), 10), degree=1)
wk_mesh

## define the spatial mesh
# TODO: currently, a fairly crude mesh is being used to decrease compute time, 
# which will bias the inferred field. at some point we should decrease the edge
# lengths to make a finer mesh
bnd <- st_as_sfc(st_bbox(loc), crs=st_crs(loc))
# ^ creates a list (?) of all points within study area, not just sampled points
loc_hex <- fm_hexagon_lattice(bnd, edge_len=0.07)
# ^ creates a hexagonal lattice of the points

sp_mesh <- fm_mesh_2d(
  loc_hex, boundary=bnd, 
  max.edge=c(0.1, 0.1), 
  crs=fm_crs(loc)
)

class(sp_mesh)
# ^ converts to a mesh which makes it an actual map (?)


# visually check the mesh looks reasonable, i.e.
# 1) relatively "even" triangles
# 2) decent amount of space around the boundries of the observation points
# 3) triangles are "small enough"
ggplot() +
  geom_fm(data=sp_mesh) +
  geom_sf(data=loc)



## define the spacetime model object
# model 121 aligns to the "critical diffusion" model (I think), a non-seperable process
st_model <- stModel.define(
  sp_mesh, wk_mesh, "121",
  control.priors = list(
    # set priors, if it matters
    prs = c(0.1, 0.05),
    prt = c(1, 0.05),
    psigma = c(1, 0.05)),
  constr = TRUE)



### Bioclim functions!!


get_raster_stack <- function(var, dates) {
  
  file_dates <- str_replace_all(dates, "-", "") |> str_sub(3,4) # only need year for now
  pattern <- str_c(var, "_equiv_20",str_c(file_dates, collapse="|"),sep="")
  files <- list.files("geofiles",  pattern=pattern, full.names=TRUE)
  as.list(stack(files))
}


extract_raster_points <- function(r, date, points) {
  rval1 <- terra::extract(r, points)
  rval1$ID <- NULL
  rvals <- as.double(unlist(rval1))
  points |> 
    st_as_sf() |> 
    mutate(rvals=rvals, year=date)
}



load_prism_var <- function(var, dates, points) {
  stack <- get_raster_stack(var, dates)
  varname <- str_replace_all(var, "-", "_")
  
  furrr::future_map2_dfr(
    stack, dates, 
    extract_raster_points, points, 
    .options=furrr_options(seed=NULL)
  ) |> 
    rename({{varname}} := rvals)
}



### Putting it all together: bioclim

biovars <- c("bio05","bio06","bio08","bio10","bio11","bio13","bio14","bio18")
bioyears <- c("2021","2022","2023","2024")

for (bvar in biovars){
  frame_bvar <- st_sf(
    id = integer(0),
    name = character(0),
    value = numeric(0),
    geometry = st_sfc(), 
    crs = st_crs(loc) # Optional: Set your Coordinate Reference System
  )
for (byear in bioyears){
  oneyear <- load_prism_var(bvar,byear,loc)
  oneyear$year <- as.double(oneyear$year)
  frame_bvar <- rbind(frame_bvar,oneyear)
} 
  traps_sf<- traps_sf %>%
    left_join(st_drop_geometry(frame_bvar),by = c("lat","long","year"))}


write.csv(traps_sf, file = "traps_sf_bc_dm.csv", row.names = TRUE)



## DAymet
tri_idx1 <- sp_mesh$graph$tv
coords1 <- sp_mesh$loc

poly_list <- lapply(1:nrow(tri_idx1), function(i) {
  # Extract the 3 vertices defining the triangle
  v_idx <- tri_idx1[i, ]
  # Close the polygon by repeating the first vertex at the end
  triangle_coords <- coords1[c(v_idx, v_idx[1]), 1:2] 
  st_polygon(list(triangle_coords))
})

mesh_spsf <- st_sf(geometry = st_sfc(poly_list, crs = crs(loc)))

st_write(mesh_spsf, "mesh_output.shp", overwrite = TRUE)

meshshape <- vect("mesh_output.shp")

srad_2024 <- rast("Daymet/Daymet_Daily/daymet_v4_daily_na_tmax_2021.nc")

meshshape <- project(meshshape,crs(srad_2024))

prcp_24 <- crop(srad_2024,meshshape)
plot(prcp_24[[140]])

writeRaster(prcp_24, "Daymet/Daymet_Daily/tmax_21.tif", overwrite = TRUE)



##Daymet helper functions

get_raster_stack2 <- function(var, dates) {
  
  pattern <- str_c("Daymet/Daymet_Daily/",var,"_",dates,".tif")
  files <- rast(pattern)
}

extract_raster_points2 <- function(r, points) {
  rval1 <- project(r,crs(points))
  time(rval1) <- as.Date(time(rval1))
  names(rval1) <-time(rval1)
  rval2 <- terra::extract(rval1, points)
  rvals <- bind_cols(st_drop_geometry(points),rval2)
  rval3 <- st_as_sf(rvals,coords = c("long","lat"),crs = crs(points))
  rval3 <- rval3 %>%
    mutate(
      long = st_coordinates(.)[, 1],
      lat = st_coordinates(.)[, 2]
    )
}

load_var2 <- function(var,dates,points){
  file <- get_raster_stack2(var,dates)
  outfile <- extract_raster_points2(file,points)
}



##,"snowweq","swradiation","vapePres" "tmax","tmin","dayl","prcp"

#Daymet add the data!!



dayvars<- c("snowweq","swradiation","vapePres")
dayyrs <- c(21,22,23,24)

for (var in dayvars){
  frame_day <- st_sf(
    lat = numeric(),
    long = numeric(),
    geometry = st_sfc(),
    var = numeric(),
    date = Date(),
    crs = st_crs(loc) # Optional: Set your Coordinate Reference System
  )
  for (year in dayyrs){
    file <- load_var2(var,year,loc)
    daydates <- as.list(colnames(file))
    daydates <- daydates[-1]
    daydates <- head(daydates,-3)
    
    
    for (day in daydates){
      df <- file %>% select("lat","long")
      df[[var]] <- file[[day]]
      df$date <- as.Date(day)
      frame_day <- rbind(frame_day,df)
      
    }
    print(year)
  }
  traps_sf<- traps_sf %>%
    left_join(st_drop_geometry(frame_day),by = c("lat","long","date"))   
   print(var)
}



traps_vars <- read.csv("traps_sf_bc_dm.csv")

## adding land cover data
lc_grid <- rast("geofiles/Annual_NLCD_LndCov_2024_CU_C1V1_mffklurgpwkbfm.tiff")

# relabel to readable land cover names, dropping any that aren't anywhere in 
# Dallas county:
values(lc_grid) <- land_cover_labels(values(lc_grid), drop=TRUE)
lc_grid <- project(lc_grid, traps_sf)

plot(lc_grid)

# model formula: just a global intercept and the field
model <- result ~ Intercept(1, prec.linear=0.5) + traps_sf$bio06 + 
   land_cover(lc_grid, model="iid", hyper=list(theta=list(prior="pc.prec", param=c(1, 0.05)))) +
  st_field(list(space=geometry, time=week_ind), model=st_model)

# fit the model
fit <- bru(model, traps_sf, family="binomial")

summary(fit)
plot(fit, "land_cover")

# make weekly prediction maps for selected time period
pred_month <- 7
pred_year <- 2024

pix_pred <- fm_pixels(sp_mesh, dims=c(100, 100), mask=bnd)

wk_pred <- filter(traps_sf, month(date) == pred_month, year == pred_year) |> 
  pull(week_ind) |> 
  unique()

pred_sf <- cross_join(pix_pred, tibble(week_ind=wk_pred))

pred_formula <- ~inla.link.invlogit(Intercept+prec_wet_mon+land_cover+st_field)
pred <- predict(fit, pred_sf, pred_formula, n.samples=1000)

ggplot() +
  gg(pred, aes(fill=mean), geom="tile") +
  geom_sf(data=bnd, fill=NA, col="black", linewidth=0.55) +
  geom_sf(aes(shape=as.factor(result)), filter(traps_sf, week_ind %in% wk_pred), col="gray70") +
  facet_wrap(~week_ind, nrow=2) +
  scale_fill_viridis_c(option="inferno") +
  scale_shape_manual(values=c(96, 4)) +
  labs(fill="prediction\nmean", shape="WNV positive") +
  theme_map()

# save as png cause it looks better with SpatialPixels
# ggsave("figs/spacetime-fit-jul-24.png", width=7.2, height=4.8, bg="white")

# examine posterior spatiotemporal hyperparameters--------------------------------
post1 <- inla.tmarginal(\(x) exp(x), fit$internal.marginals.hyperpar[[1]]) |> 
  as_tibble() |> 
  mutate(parameter="spatial range")

post2 <- inla.tmarginal(\(x) exp(x), fit$internal.marginals.hyperpar[[2]]) |> 
  as_tibble() |> 
  mutate(parameter="temporal range")

post3 <- inla.tmarginal(\(x) exp(x), fit$internal.marginals.hyperpar[[3]]) |> 
  as_tibble() |> 
  mutate(parameter="field std.")

bind_rows(post1, post2, post3) |> 
  mutate(parameter=fct_inorder(parameter)) |> 
  ggplot(aes(x, y, col=parameter)) +
  geom_line() +
  facet_wrap(~parameter, scales="free")

# confirm unvisited land cover classes are predicted correctly (Bren TODO)--------
# pred_lc <- pix_pred |> 
#     mutate(land_cover=extract(lc_grid, pix_pred)[,2])
# 
# lc_classes_pred_grid <- unique(extract(lc_grid, pix_pred)[,2])
# 
# pred_lc <- predict(fit, tibble(land_cover=lc_classes_pred_grid), ~land_cover, n.samples=500)
# 
# pred_lc
# 
# ggplot(, aes(land_cover))
# 
# pred_lc |> 
#     filter(`land_cover.ID` %in% unique(fit$summary.random$land_cover$ID))