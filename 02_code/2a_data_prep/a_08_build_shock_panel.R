############ SHOCK / BASELINE EXPOSURE PANEL ############
# Day-level version of the trend-vs-shock decomposition in Zou, Carlton & Grover (medRxiv
# 2026.03.30.26349797). They work at county-month and split each weather variable into a
# 90-day rolling baseline (sustained conditions) and a shock (recent period minus that
# baseline). Same idea here, daily.
#
# Why bother: the case-crossover has been fitting exposure LEVELS, and a level is not the
# same quantity as a departure from a cell's own recent baseline. Soil moisture 20% below
# where it's been sitting for three months is a drought signal; soil moisture at 0.25 m3/m3
# is just a number.
#
# One structural caveat, measured not assumed: the baselines are nearly useless in this
# design. Share of variance surviving within-month stratum centring:
#   runoff shock 85% | precipitation shock 73% | soil moisture shock 36%
#   runoff base90 12% | temperature shock 11% | precipitation base90 9%
#   temperature base90 2%
# Anything on a seasonal timescale gets eaten by the strata, so only the shocks are
# estimable here. Sustained conditions need the cohort model, not this one.

rm(list = ls())
project.folder = paste0(print(here::here()), '/')
source(paste0(project.folder, 'create_folder_structure.R'))
source(paste0(functions.folder, 'script_initiate.R'))

BASE_WIN <- 90L        # rolling baseline, days
RECENT   <- 7L         # "recent conditions" window the shock is measured against
LAGS     <- 0:28       # lag columns to write out, matching the met crossbases
# wind belongs here too: Zou et al. carry high-wind coverage as a storm predictor, and it was
# in the level model, so leaving it out of the shock set was an oversight rather than a choice.
SH_VARS  <- c("temperature", "precipitation", "soil_moisture", "runoff",
             "wind_speed", "relative_humidity", "snow_cover")

# 1. Daily source series, one row per cell-day across all 7 states.
daily <- setDT(readRDS(paste0(objects_folder, multistate_ts_daily_rds)))
setorder(daily, zone_id, date)
daily[, date := as.Date(date)]

# 1a. ERA5 ships runoff and precipitation in METRES (order 1e-4). Convert before any
#     accumulation or the shocks come out at a scale nothing else in the pipeline uses.
for (v in intersect(c("precipitation", "runoff"), SH_VARS))
  set(daily, j = v, value = daily[[v]] * 1000)

# 2. Baseline and shock per variable.
#    base90  = mean over the 90 days ENDING YESTERDAY (shift(.,1) keeps day t out of its own
#              baseline, same temporal-leakage guard Zou et al. use)
#    recent7 = mean over the last 7 days including today
#    shock   = recent7 - base90
# 2a. Precipitation and runoff get log1p FIRST, then differenced, so the shock is a
#     proportional departure rather than a difference in mm. You can't log a shock after the
#     fact (it goes negative), and the raw differences are badly skewed: precipitation
#     +1.00 -> +0.25 and runoff +4.65 -> +2.07 once built on the log scale. Temperature and
#     soil moisture are near-symmetric already (-0.37, -0.64), so leave them alone.
LOG_SHOCK <- c("precipitation", "runoff")
for (v in SH_VARS) {
  x <- if (v %in% LOG_SHOCK) log1p(pmax(daily[[v]], 0)) else daily[[v]]
  set(daily, j = "..x", value = x)
  daily[, (paste0(v, "_base90"))  := shift(frollmean(..x, BASE_WIN, align = "right"), 1),
        by = zone_id]
  daily[, (paste0(v, "_recent7")) := frollmean(..x, RECENT, align = "right"), by = zone_id]
  daily[, (paste0(v, "_shock"))   := get(paste0(v, "_recent7")) - get(paste0(v, "_base90"))]
}
daily[, ..x := NULL]

SHOCK <- paste0(SH_VARS, "_shock")
ex <- daily[, c("zone_id", "date", SHOCK), with = FALSE]
saveRDS(ex, paste0(objects_folder, "daily_shock_predictors.RDS"))
cat("first date with a full 90-day baseline:",
    as.character(min(ex[!is.na(temperature_shock)]$date)), "\n")

# 3. Attach shock_Lag0..Lag28 to the case-crossover panel.
cc <- setDT(readRDS(paste0(objects_folder, "case_crossover_df_timestrat_month_post7.RDS")))
# a_08 emits `date` twice (identical copies). harmless until you try to merge, then fatal.
if (anyDuplicated(names(cc))) cc <- cc[, which(!duplicated(names(cc))), with = FALSE]
cc[, date := as.Date(date)]

# one lag at a time: joining date - L back onto the daily series never materialises 29 copies
# of a 4.9M-row table, which is what happens if you shift in place
for (L in LAGS) {
  src <- ex[, c(.(zone_id = zone_id, date = date + L), .SD), .SDcols = SHOCK]
  setnames(src, SHOCK, paste0(SHOCK, "_Lag", L))
  cc <- merge(cc, src, by = c("zone_id", "date"), all.x = TRUE, sort = FALSE)
}

# 4. z-score each lag column, matching how a_08 scaled the raw exposures.
for (v in SHOCK) for (L in LAGS) {
  cn <- paste0(v, "_Lag", L)
  set(cc, j = cn, value = as.numeric(scale(cc[[cn]])))
}

# 5. Drop rows without a complete lag history, then drop strata the deletions broke. A
#    stratum with no case, or a case with no referents, contributes nothing to the
#    conditional likelihood and will silently distort the counts if left in.
need <- paste0(rep(SHOCK, each = length(LAGS)), "_Lag", LAGS)
cat(sprintf("rows %d -> %d complete | cases %d -> %d\n", nrow(cc), sum(complete.cases(cc[, ..need])),
            sum(cc$outbreak_binary), sum(cc$outbreak_binary[complete.cases(cc[, ..need])])))
cc <- cc[complete.cases(cc[, ..need])]
cc[, `:=`(nc = sum(outbreak_binary), nr = sum(outbreak_binary == 0)), by = stratum_id]
cc <- cc[nc == 1L & nr > 0L][, c("nc", "nr") := NULL]

saveRDS(cc, paste0(objects_folder, "case_crossover_df_shock_post7.RDS"))
cat(sprintf("shock panel: cases %d (MN %d) | refs %d | %.1f per case\n",
            sum(cc$outbreak_binary), sum(cc[state == "Minnesota"]$outbreak_binary),
            sum(cc$outbreak_binary == 0),
            sum(cc$outbreak_binary == 0) / sum(cc$outbreak_binary)))
