library(tidyverse)

land_cover_labels <- function(lc, drop=FALSE) {
    # lc_map <- as.character(c(11, 12, 21:24, 31, 41:43, 51, 52, 71:74, 81, 82, 90, 95))
    # names(lc_map) <- land_cover_names()
    # ret <- fct_recode(lc, !!!lc_map)
    ret <- recode_values(
        as.character(lc),
        from=as.character(c(11, 12, 21:24, 31, 41:43, 51, 52, 71:74, 81, 82, 90, 95)),
        to=land_cover_names()
    ) |> 
        factor()
    
    if (drop)
        return(fct_drop(ret))
    ret
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