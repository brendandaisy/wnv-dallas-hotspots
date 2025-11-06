# prep-all-years.R----------------------------------------------------------------
# prepare the data from multiple years (currently 2021-2024) into a single--------
# shapefile with lat-long found using tidygeocoder--------------------------------
# --------------------------------------------------------------------------------
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

# assumes the tables are in folder 'comprehensive files'
# TODO: ASSUMPTION: currently assumes any traps with zero mosquitoes are either 
# a false for WNV, unless either `qk_total` or `Result` is marked with an "-", in which case
# assume the trap was not visited / NA. We should check with Dallas Co about this.
qk21 <- readxl::read_xlsx("comprehensive files/Comprehensive_WNV_DYS_tables_2021.xlsx", sheet="G") |> 
    transmute(
        date=as.Date(`Collected`),
        epiweek=epiweek(date),
        zone=Zone,
        address=Address,
        city=City,
        zip=Zip,
        # country="USA",
        num_trapped=qk_total,
        result=case_match(
            Result,
            c("N", "ID") ~ FALSE,
            c("Pos", "POS") ~ TRUE,
            .default=NA
        )
    )

qk22 <- readxl::read_xlsx("comprehensive files/Comprehensive_WNV_DYS_tables_2022.xlsx", sheet="G") |>
    transmute(
        date=as.Date(`Collected`),
        epiweek=epiweek(date),
        zone=Zone,
        address=Address,
        city=City,
        zip=Zip,
        num_trapped=as.character(qk_total),
        result=case_match(
            Result,
            c("N", "ID") ~ FALSE,
            c("Pos", "POS") ~ TRUE,
            .default=NA
        )
    )

qk23 <- readxl::read_xlsx("comprehensive files/Comprehensive_WNV_DYS_tables_2023.xlsx", sheet="G") |>
    transmute(
        date=as.Date(`Collected`),
        epiweek=epiweek(date),
        zone=Zone,
        address=Address,
        city=City,
        zip=as.character(Zip),
        # country="USA",
        num_trapped=qk,
        result=case_match(
            Result,
            c("N", "ID") ~ FALSE,
            c("Pos", "POS") ~ TRUE,
            .default=NA
        )
    )

qk24 <- readxl::read_xlsx("comprehensive files/Comprehensive_WNV_DYS_tables_2024.xlsx", sheet="G") |>
    transmute(
        date=as.Date(`Collected`),
        epiweek=epiweek(date),
        zone=Zone,
        address=Address,
        city=City,
        zip=as.character(Zip),
        num_trapped=qk,
        result=case_match(
            Result,
            c("N", "ID") ~ FALSE,
            c("Pos", "POS") ~ TRUE,
            .default=NA
        )
    )

# combine the dataframes and filter any final NAs
traps0 <- bind_rows(qk21, qk22, qk23, qk24) |> 
    mutate(
        country="USA", 
        num_trapped=parse_number(str_replace(num_trapped, "-", "0"))
    ) |> 
    filter(!is.na(result), !is.na(num_trapped))

# found that 'arcgis' had the most success and was able to find all coords in one go:
traps_geo <- geocode(traps0, method="arcgis", street=address, city=city, postalcode=zip, country=country)

# convert to 'sf' format
traps_sf <- traps_geo |> 
    st_as_sf(coords=c("long", "lat"), crs="WGS84", remove=FALSE)

# write as a shapefile
write_sf(traps_sf, "traps-all-years.shp")

## Some visualizations-------------------------------------------------------------
library(zoo)

traps_date <- traps_sf |> 
    st_drop_geometry() |> 
    group_by(date) |> 
    summarize(num_pos=sum(result)) |> 
    mutate(pos_wk=rollapply(num_pos, 7, mean, fill=NA))

ggplot(traps_date, aes(date)) +
    geom_point(aes(y=num_pos), alpha=0.2) +
    geom_line(aes(y=pos_wk), col="blue") +
    theme_half_open()

## number of visits each week
traps_sf |> 
    st_drop_geometry() |> 
    count(epiweek, year=year(date)) |> 
    ggplot(aes(interaction(epiweek, year), n)) +
    geom_col() +
    theme_half_open()
