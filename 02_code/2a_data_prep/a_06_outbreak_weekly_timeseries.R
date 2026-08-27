## ============================================================
## Build a long outbreak × week panel covering 52 ISO weeks of 2022
## plus any prior weeks (late Nov / Dec 2021) required as
## case-crossover controls for early-January 2022 cases.
##
## Output columns:
##   outbreak_id        — source CSV identifier
##   premises           — farm identifier (same source CSV)
##   state              — 2-letter state code
##   state_name         — full state name (e.g. "Minnesota")
##   flyway             — ATL / CEN / MISS / PAC (from source CSV)
##   herd_flock_type    — Commercial / Backyard
##   production_type    — finer type (commercial turkey, etc.)
##   link               — IND / INR / CSLT / NOSEQ (for downstream filtering)
##   lon, lat           — outbreak coordinates (degrees, WGS84)
##   week_id            — "YYYY-Www" ISO week tag
##   iso_year, iso_week — separated for ordering
##   week_start         — Monday of that ISO week (Date)
##   case_week_start    — Monday of the outbreak's own case-week (Date)
##   is_case_week       — 1 if this row is the outbreak's case-week, 0 otherwise
##   lag_weeks_from_case— signed integer (case-week = 0; control candidates are -1..-4)
##
## Save targets:
##   01_data/1b_outcome_data/national_spillover_events_casecrossover_2022.csv
##   01_data/1b_outcome_data/national_spillover_events_casecrossover_2022.gpkg
## ============================================================

rm(list = ls())

# 0a. Declare root directory, folder locations and load essential stuff
project.folder <- paste0(print(here::here()), "/")
source(paste0(project.folder, "create_folder_structure.R"))
source(paste0(functions.folder, "script_initiate.R"))

suppressPackageStartupMessages({
  library(data.table); library(lubridate); library(sf)
})

# 1. Load source outbreak data
src <- fread(paste0(spillover_data, "hpai-outbreaks-clean-placeholder-032124.csv"))

# 2. Restrict to 2022 outbreaks. Keep ALL link categories so downstream code
# can filter (IND/INR for primary, all for sensitivity comparisons).
src[, date_confirmed := as.Date(date_confirmed)]
src[, year := lubridate::year(date_confirmed)]
out22 <- src[year == 2022 & !is.na(lon) & !is.na(lat)]

# Add full state name via base R's state.abb / state.name lookup
state_lookup <- data.table(state = c(state.abb, "DC"),
                            state_name = c(state.name, "District of Columbia"))
out22 <- merge(out22, state_lookup, by = "state", all.x = TRUE)

# 3. Compute each outbreak's case-week (ISO week, Monday start)
out22[, case_week_start := lubridate::floor_date(date_confirmed, unit = "week",
                                                  week_start = 1)]
out22[, iso_year_case   := lubridate::isoyear(date_confirmed)]
out22[, iso_week_case   := lubridate::isoweek(date_confirmed)]

# 4. Build the shared weekly calendar
# The panel needs to support both the WEEKLY case-crossover (case-week W vs
# controls at W-1..W-4) and the DAILY case-crossover (each day in that
# window plus a daily lag-exposure history). For the earliest 2022 case
# the back-look is:
#       4 weeks  (case-crossover control window: W-1..W-4)
#     + 8 weeks  (daily lag-56 buffer; covers lag-28 / -42 / -56 sensitivities)
#     = 12 ISO weeks before the earliest case-week start.
#
# Panel always covers ISO weeks 1..52 of 2022 at minimum. The backward
# extension into Nov–Dec 2021 lets eBird and meteorology data join cleanly
# onto every required day of the daily case-crossover, including the longer
# lag-window sensitivity analyses.
iso_week1_2022  <- as.Date("2022-01-03")          # Monday of ISO week 1, 2022
iso_week52_2022 <- as.Date("2022-12-26")          # Monday of ISO week 52, 2022
daily_lag_buffer_weeks   <- 8                     # 56 d ÷ 7 d/wk = 8 weeks
crossover_control_weeks  <- 4                     # W-1..W-4 controls
backlook_weeks <- crossover_control_weeks + daily_lag_buffer_weeks  # = 12
panel_start <- min(iso_week1_2022,
                   min(out22$case_week_start) -
                     lubridate::weeks(backlook_weeks))
panel_end   <- iso_week52_2022
weeks_grid <- data.table(
  week_start = seq(panel_start, panel_end, by = "1 week")
)
weeks_grid[, `:=`(
  iso_year = lubridate::isoyear(week_start),
  iso_week = lubridate::isoweek(week_start),
  week_id  = sprintf("%d-W%02d",
                     lubridate::isoyear(week_start),
                     lubridate::isoweek(week_start))
)]

cat("Panel covers ", nrow(weeks_grid), " weeks ",
    "(", as.character(panel_start), " to ", as.character(panel_end), ")\n", sep = "")
cat("2022 outbreaks in panel: ", nrow(out22), "\n", sep = "")

# 5. Cross-join outbreaks × weeks (long panel)
ts_panel <- CJ(outbreak_id = out22$outbreak_id, week_start = weeks_grid$week_start,
               sorted = FALSE)
# Pull outbreak-level columns onto the panel
out22_cols <- c("outbreak_id", "premises", "state", "state_name", "flyway",
                "herd_flock_type", "production_type", "link",
                "lon", "lat", "case_week_start")
ts_panel <- merge(ts_panel, out22[, ..out22_cols], by = "outbreak_id", all.x = TRUE)
# Pull week-level columns
ts_panel <- merge(ts_panel,
                  weeks_grid[, .(week_start, iso_year, iso_week, week_id)],
                  by = "week_start", all.x = TRUE)

# 6. Case-week flag and lag-from-case (in weeks)
ts_panel[, is_case_week        := as.integer(week_start == case_week_start)]
ts_panel[, lag_weeks_from_case := as.integer((week_start - case_week_start) / 7)]

# 7. Order columns for readability and sort
setcolorder(ts_panel, c(
  "outbreak_id", "premises", "state", "state_name", "flyway",
  "herd_flock_type", "production_type", "link",
  "lon", "lat",
  "week_id", "iso_year", "iso_week", "week_start",
  "case_week_start", "is_case_week", "lag_weeks_from_case"
))
setorder(ts_panel, outbreak_id, week_start)

# 8a. Save the long panel as CSV
out_csv <- paste0(spillover_data, "national_spillover_events_casecrossover_2022.csv")
fwrite(ts_panel, out_csv)

# 8b. Save the same panel as a GPKG (WGS84 point geometry from lon/lat).
# Dates need to be coerced to character for GPKG compatibility (the OGR
# Date driver round-trips fine but mixed Date columns sometimes confuse
# downstream readers — store as ISO strings).
ts_sf <- st_as_sf(ts_panel, coords = c("lon", "lat"), crs = 4326, remove = FALSE)
# Coerce Date columns to character (ISO yyyy-mm-dd) for portability
for (col in c("week_start", "case_week_start")) {
  if (col %in% names(ts_sf)) ts_sf[[col]] <- as.character(ts_sf[[col]])
}
out_gpkg <- paste0(spillover_data, "national_spillover_events_casecrossover_2022.gpkg")
if (file.exists(out_gpkg)) file.remove(out_gpkg)
st_write(ts_sf, out_gpkg,
         layer = "national_spillover_events_casecrossover_2022",
         quiet = TRUE)

cat("\nSaved:\n",
    "  CSV : ", out_csv, "\n",
    "  GPKG: ", out_gpkg, "\n", sep = "")
cat("Rows : ", nrow(ts_panel),
    " (", nrow(out22), " outbreaks × ", nrow(weeks_grid), " weeks)\n", sep = "")
cat("Case-week rows (one per outbreak): ", sum(ts_panel$is_case_week), "\n", sep = "")
cat("Control-window rows (lag -1..-4 weeks): ",
    sum(ts_panel$lag_weeks_from_case %in% -4:-1), "\n", sep = "")
cat("Flyways present: ",
    paste(sort(unique(ts_panel$flyway)), collapse = ", "), "\n", sep = "")
