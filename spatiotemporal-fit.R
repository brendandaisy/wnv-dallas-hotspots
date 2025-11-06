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

traps_sf <- read_sf("traps-all-years.shp") |> 
    rename(num_trapped=nm_trpp) |> 
    mutate(year=year(date), .after=epiweek) |> 
    group_by(year, epiweek) |> 
    mutate(week_ind=cur_group_id())

loc <- distinct(traps_sf, lat, long, geometry)
wk <- seq(min(traps_sf$week_ind), max(traps_sf$week_ind), 1)

## define the temporal mesh
wk_mesh <- fm_mesh_1d(wk, degree=1)
wk_mesh

## define the spatial mesh
# TODO: currently, a fairly crude mesh is being used to decrease compute time, 
# which will bias the inferred field. at some point we should decrease the edge
# lengths to make a finer mesh
bnd <- st_as_sfc(st_bbox(loc), crs=st_crs(loc))
loc_hex <- fm_hexagon_lattice(bnd, edge_len=0.05)

sp_mesh <- fm_mesh_2d(
    loc_hex, boundary=bnd, 
    max.edge=c(0.07, 0.1), 
    crs=fm_crs(loc)
)

ggplot() +
    geom_fm(data=sp_mesh) +
    geom_sf(data=loc)

# define the spacetime model object
# model 121 aligns to the "critical diffusion" model (I think), a non-seperable process
stm121 <- stModel.define(
    sp_mesh, wk_mesh, "121",
    control.priors = list(
        # set priors, if it matters
        prs = c(0.2, 0.05),
        prt = c(5, 0.05),
        psigma = c(1, 0.05)),
    constr = TRUE)

# model formula: just a global intercept and the field
model <- result ~ Intercept(1) + 
    field(list(space=geometry, time=week_ind), model=stm121)

# fit the model
fit <- bru(model, traps_sf, family="binomial")

summary(fit)

# make weekly predictions for a dense grid of points from 2024
pix_pred <- fm_pixels(sp_mesh, dims=c(150, 150), mask=bnd)
wk_pred <- filter(traps_sf, month(date) == 7, year == 2024)$week_ind |> unique()
pred_sf <- cross_join(pix_pred, tibble(week_ind=wk_pred))
    
pred <- predict(fit, pred_sf, ~inla.link.invlogit(field + Intercept), n.samples=1000)

ggplot() +
    gg(pred, aes(fill=mean), geom="tile") +
    geom_sf(data=bnd, fill=NA, col="black", linewidth=0.55) +
    geom_sf(aes(shape=as.factor(result)), filter(traps_sf, week_ind %in% wk_pred), col="gray70") +
    facet_wrap(~week_ind, nrow=2) +
    scale_fill_viridis_c(option="inferno") +
    scale_shape_manual(values=c(1, 4)) +
    labs(fill="prediction\nmean", shape="WNV positive") +
    theme_map()

# save as png cause it looks better with SpatialPixels
ggsave("figs/spacetime-fit-jul-24.png", width=7.2, height=4.8, bg="white")

#  examine posterior hyperparameters----------------------------------------------
post1 <- inla.tmarginal(\(x) exp(x), fit$internal.marginals.hyperpar[[1]]) |> 
    as_tibble() |> 
    mutate(parameter="spatial range")

# TODO: idk why tf the temporal range is so huge, but its causing bad predictions
# during non-summer months, plus just makes no sense
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
