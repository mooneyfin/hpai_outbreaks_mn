# a_00_era5_land_rgee.R
# R-native ERA5-Land daily extraction via the rgee package.
# Mirrors a_00_era5_land_gee.js (Google Earth Engine JavaScript editor) so the
# JS and R versions are interchangeable — pick whichever you prefer to debug.
#
# Why have both: the JS editor in code.earthengine.google.com is the easiest
# place to inspect intermediate results visually; the R version lets the daily
# extraction live inside the SPARK Lab project structure and pull straight into
# tibbles without the Drive → download → CSV step.
#
# Requirements:
#   - install.packages('rgee'); rgee::ee_install() one-time Python setup
#   - rgee::ee_Initialize(drive = TRUE) one-time GEE auth (uses your Google account)
#   - The fishnet asset projects/ee-feedlot/assets/mn_fishnet_10km_final must exist
#     (created from the bottom of a_00_era5_land_gee.js the first time)

# 0a. Declare root directory, folder locations, and load packages + helpers
rm(list = ls())
project.folder = paste0(print(here::here()), '/')
source(paste0(project.folder, 'create_folder_structure.R'))
source(paste0(functions.folder, 'script_initiate.R'))

# 0b. rgee initialization (interactive — first run will open browser auth)
library(rgee)
ee_Initialize(drive = TRUE)

# 1a. Load assets
fishnet  <- ee$FeatureCollection('projects/ee-feedlot/assets/mn_fishnet_10km_final')
era5land <- ee$ImageCollection("ECMWF/ERA5_LAND/DAILY_AGGR")

# 1b. RH and wind helpers — same expressions as the JS version
compute_rh <- function(t_k, td_k) {
  t_c  <- t_k$subtract(273.15)
  td_c <- td_k$subtract(273.15)
  es <- t_c$expression(
    '0.6108 * exp((17.27 * T) / (T + 237.3))', list(T = t_c)
  )
  ea <- td_c$expression(
    '0.6108 * exp((17.27 * Td) / (Td + 237.3))', list(Td = td_c)
  )
  ea$divide(es)$multiply(100)$rename('relative_humidity')
}
compute_wind <- function(u, v) {
  u$pow(2)$add(v$pow(2))$sqrt()$rename('wind_speed')
}

# 1c. Build a single-day zonal feature collection
daily_zonal <- function(date_str) {
  day      <- ee$Date(date_str)
  next_day <- day$advance(1, 'day')
  era_day  <- era5land$filterDate(day, next_day)

  # DAILY_AGGR already aggregates over the day; .mean() of one image returns
  # that image's values unchanged. Kept for parity with the JS version.
  mean_vars <- era_day$mean()
  sum_vars  <- era_day$select(c(
    'total_precipitation_sum', 'runoff_sum', 'total_evaporation_sum'
  ))$mean()

  t_k  <- mean_vars$select('temperature_2m')
  td_k <- mean_vars$select('dewpoint_temperature_2m')
  rh   <- compute_rh(t_k, td_k)
  wind <- compute_wind(
    mean_vars$select('u_component_of_wind_10m'),
    mean_vars$select('v_component_of_wind_10m')
  )

  combined <- t_k$rename('temperature')$
    addBands(td_k$rename('dewpoint'))$
    addBands(rh)$
    addBands(wind)$
    addBands(mean_vars$select('volumetric_soil_water_layer_1')$rename('soil_moisture'))$
    addBands(mean_vars$select('snow_density'))$
    addBands(sum_vars$select('total_precipitation_sum')$rename('precipitation'))$
    addBands(sum_vars$select('runoff_sum'))$
    addBands(sum_vars$select('total_evaporation_sum')$rename('evapotranspiration'))$
    addBands(mean_vars$select('leaf_area_index_high_vegetation')$rename('leaf_area_index'))$
    addBands(mean_vars$select('surface_net_solar_radiation_sum')$rename('shortwave_radiation'))$
    addBands(mean_vars$select('surface_thermal_radiation_downwards_sum')$rename('longwave_radiation'))$
    addBands(mean_vars$select('surface_pressure'))

  # Pull zonal means per fishnet cell
  stats <- combined$reduceRegions(
    collection = fishnet,
    reducer    = ee$Reducer$mean(),
    scale      = 10000
  )$map(function(f) f$set('date', date_str))

  stats
}

# 2a. Date range — Nov 15 2021 through Dec 31 2022 (lag-28 windows of early-2022
#     cases reach into Dec 2021; see code comments in a_03_create_casecrossover_df.Rmd)
date_seq <- format(
  seq.Date(as.Date("2021-11-15"), as.Date("2022-12-31"), by = "day"),
  "%Y-%m-%d"
)

# 2b. Extract one day at a time and bind into an R tibble. For the full
#     ~412-day range × ~5500 zones this fetches ~2.3 M rows; doable in one
#     session over a residential connection. For batch runs, use the
#     Drive-export branch below instead.
cat("Fetching", length(date_seq), "days of zonal ERA5-Land...\n")
daily_list <- vector("list", length(date_seq))
for (i in seq_along(date_seq)) {
  d <- date_seq[[i]]
  daily_list[[i]] <- ee_as_sf(daily_zonal(d), maxFeatures = 10000) |>
    sf::st_drop_geometry() |>
    tibble::as_tibble()
  if (i %% 30 == 0) cat("  ", i, "/", length(date_seq), "days\n")
}
daily_df <- dplyr::bind_rows(daily_list)

# 3a. Clean up snow_density artefact (ERA5-Land returns ~99.99998 when no snow)
daily_df <- daily_df |>
  dplyr::mutate(snow_density = ifelse(snow_density <= 100, 0, snow_density))

# 3b. Save to the exposure-data folder
saveRDS(daily_df,
        file.path(env_data, "mn_daily_era5land_2021_2022.rds"))
data.table::fwrite(daily_df,
                   file.path(env_data, "mn_daily_era5land_2021_2022.csv"))
cat("Wrote", nrow(daily_df), "rows to", env_data, "\n")


# ─────────────────────────────────────────────────────────────────────────────
# OPTIONAL: Drive-export branch for unattended batch runs
# ─────────────────────────────────────────────────────────────────────────────
# For very large fetches (e.g., 1997–2022 weekly climatology baselines), the
# in-session pull above will be too slow. Use ee_table_to_drive() to kick off
# async exports the same way the JavaScript script does, then download the
# CSVs from Google Drive once they finish.
#
# Example for one year:
#
#   year_fc <- ee$FeatureCollection(lapply(
#     format(seq.Date(as.Date("2022-01-01"), as.Date("2022-12-31"), by = "day"),
#            "%Y-%m-%d"),
#     daily_zonal
#   ))$flatten()
#
#   task <- ee_table_to_drive(
#     collection     = year_fc,
#     description    = "MN_ERA5Land_Daily_2022_rgee",
#     folder         = "EarthEngine_MN",
#     fileNamePrefix = "mn_daily_era5land_2022_rgee",
#     fileFormat     = "CSV"
#   )
#   task$start()
#   ee_monitoring(task)
