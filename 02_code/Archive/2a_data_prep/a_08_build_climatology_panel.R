############ DAILY CLIMATOLOGICAL ANOMALY PANEL ############
# A third exposure parameterisation, all seven flyway states, DAILY resolution:
#
#     anomaly = (daily 2022 value - 25-year mean for that week) / 25-year sd for that week
#
# The climatology is weekly because that is what the GEE export produces (52 week_idx values
# per cell), but the numerator stays daily. So a cell-day is compared against what that week
# of the year normally looks like in that cell, which is a proper standardised anomaly and
# keeps the daily lag structure the rest of the pipeline uses.
#
# Why this is not the 90-day rolling anomaly relabelled: a 90-day baseline barely moves
# inside a one-month stratum, so the stratum intercept absorbs it and the "anomaly" ends up
# correlating 0.87-0.99 with the raw level. A climatological normal changes week to week, so
# it survives the centring - within-stratum correlation with the level is only 0.19-0.40.

rm(list = ls())
project.folder = paste0(print(here::here()), '/')
source(paste0(project.folder, 'create_folder_structure.R'))
source(paste0(functions.folder, 'script_initiate.R'))

VARS <- c("temperature", "soil_moisture", "runoff", "wind_speed", "precipitation")
LAGS <- 0:28

# 1. Climatology: four chunks of 13 weeks each, 13,920 cells x 52 weeks.
need <- c("zone_id", "week_idx", paste0(VARS, "_mean"), paste0(VARS, "_sd"))
clim <- rbindlist(lapply(
  list.files(env_data, pattern = "flyway8_climatology_1997_2021_w", full.names = TRUE),
  fread, select = need))
cat(sprintf("climatology: %s rows, %d cells, %d weeks\n",
            format(nrow(clim), big.mark = ","), uniqueN(clim$zone_id), uniqueN(clim$week_idx)))
setkey(clim, zone_id, week_idx)

# 2. Daily observed. week_idx is 7-day blocks from Jan 1, matching the GEE definition
#    exactly - ISO weeks would land week_idx on a different calendar week each year.
daily <- setDT(readRDS(paste0(objects_folder, multistate_ts_daily_rds)))
daily[, date := as.Date(date)]
daily[, week_idx := pmin(as.integer(strftime(date, "%j")) %/% 7L, 51L)]

# The daily panel and the climatology have to be on the same units before differencing:
# ERA5 ships depths in metres in both, so nothing to convert, but temperature is Kelvin in
# both too. Leave everything as exported.
d <- merge(daily[, c("zone_id", "date", "week_idx", "state", VARS), with = FALSE],
           clim, by = c("zone_id", "week_idx"), all.x = TRUE)

# 3. z = (observed - weekly mean) / weekly sd. A zero sd means the variable never varies in
#    that week across 25 years (snow in July, say); those become NA rather than Inf.
for (v in VARS) {
  m <- d[[paste0(v, "_mean")]]; s <- d[[paste0(v, "_sd")]]
  set(d, j = paste0(v, "_clim"), value = fifelse(s > 0, (d[[v]] - m) / s, NA_real_))
}
ex <- d[, c("zone_id", "date", paste0(VARS, "_clim")), with = FALSE]
saveRDS(ex, paste0(objects_folder, "daily_climatology_anomaly.RDS"))
cat(sprintf("anomalies built for %s cell-days (%.1f%% non-missing)\n",
            format(nrow(ex), big.mark = ","),
            100 * mean(!is.na(ex$temperature_clim))))

# 4. Attach lag columns to the case-crossover panel, same way a_07 does it
CL <- paste0(VARS, "_clim")
cc <- setDT(readRDS(paste0(objects_folder, "case_crossover_df_timestrat_month_post7.RDS")))
if (anyDuplicated(names(cc))) cc <- cc[, which(!duplicated(names(cc))), with = FALSE]
cc[, date := as.Date(date)]
for (L in LAGS) {
  src <- ex[, c(.(zone_id = zone_id, date = date + L), .SD), .SDcols = CL]
  setnames(src, CL, paste0(CL, "_Lag", L))
  cc <- merge(cc, src, by = c("zone_id", "date"), all.x = TRUE, sort = FALSE)
}
for (v in CL) for (L in LAGS) {
  cn <- paste0(v, "_Lag", L); set(cc, j = cn, value = as.numeric(scale(cc[[cn]])))
}
keep <- paste0(rep(CL, each = length(LAGS)), "_Lag", LAGS)
cat(sprintf("rows %d -> %d complete | cases %d -> %d\n", nrow(cc),
            sum(complete.cases(cc[, ..keep])), sum(cc$outbreak_binary),
            sum(cc$outbreak_binary[complete.cases(cc[, ..keep])])))
cc <- cc[complete.cases(cc[, ..keep])]
cc[, `:=`(nc = sum(outbreak_binary), nr = sum(outbreak_binary == 0)), by = stratum_id]
cc <- cc[nc == 1L & nr > 0L][, c("nc", "nr") := NULL]

saveRDS(cc, paste0(objects_folder, "case_crossover_df_climanom_post7.RDS"))
cat(sprintf("climatology panel: %d cases (MN %d) | %.1f referents per case\n",
            sum(cc$outbreak_binary), sum(cc[state == "Minnesota"]$outbreak_binary),
            sum(cc$outbreak_binary == 0) / sum(cc$outbreak_binary)))
