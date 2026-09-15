############ 25-YEAR CLIMATOLOGICAL ANOMALY (MINNESOTA, WEEKLY) ############
# A third exposure parameterisation: (2022 weekly value - 1997-2021 weekly mean) / 25-yr sd.
#
# Why it is not just the level model relabelled. A 90-day rolling baseline barely moves
# inside a one-month stratum, so the stratum intercept absorbs it and the "anomaly" ends up
# correlating 0.87-0.99 with the raw level. A climatological normal moves week to week, so it
# survives the centring: within-stratum correlation with the level is only 0.19-0.40.
#
# Coverage is the constraint. The 25-year climatology was exported for Minnesota only
# (a_00_era5_land_gee_mn8_climatology.js sets TARGET_STATE = 'Minnesota'), and it is weekly.
# So this is an MN-only, weekly sensitivity - no pooled panel. That is on reasonable footing
# because the weekly model already reproduces the daily findings (Table S6).

rm(list = ls())
project.folder = paste0(print(here::here()), '/')
source(paste0(project.folder, 'create_folder_structure.R'))
source(paste0(functions.folder, 'script_initiate.R'))
suppressMessages({library(INLA); library(Matrix)})
INLA::inla.setOption(num.threads = 1)
source(paste0(functions.folder, 'inla_dlnm_helpers.R'))

PREC <- 4; NS2 <- list(fun = "ns", df = 2); LIN <- list(fun = "lin")
SHV  <- c("temperature", "soil_moisture", "runoff", "wind_speed")
LABS <- c(temperature = "Temperature", soil_moisture = "Soil moisture",
          runoff = "Runoff", wind_speed = "Wind speed", anseriformes = "Anseriformes")
cell <- function(e) sprintf("%.2f (%.2f, %.2f)", e[["OR"]], e[["low"]], e[["high"]])

# 1. Panel. Use the shock panel so all three parameterisations sit on identical rows.
D <- setDT(readRDS(paste0(objects_folder, "case_crossover_df_shock_post7.RDS")))
if (anyDuplicated(names(D))) D <- D[, which(!duplicated(names(D))), with = FALSE]
D <- D[state == "Minnesota"]
D[, date := as.Date(date)]

# 2. Attach the weekly z-anomaly. The climatology is keyed on week_start, so match each
#    cell-day to the week containing it. The panel starts 2021-11-15 but the climatology
#    covers 2022 only, so pre-2022 rows have to be checked rather than assumed.
an <- setDT(readRDS(paste0(objects_folder, "climate_anomaly_panel.RDS")))
ZC <- paste0(SHV, "_z_25yr")
an <- an[, c("zone_id", "week_start", ZC), with = FALSE]
an[, week_start := as.Date(week_start)]
setkey(an, zone_id, week_start)

lookup <- function(zid, dts) {
  # roll = TRUE carries the most recent preceding week forward; rows before the first
  # climatology week get NA rather than a silently wrong value
  m <- an[.(zid, dts), roll = TRUE, on = .(zone_id, week_start)]
  m[dts < min(an$week_start), (ZC) := NA]
  m
}
for (L in 0:28) {
  m <- lookup(D$zone_id, D$date - L)
  for (v in SHV) set(D, j = paste0(v, "_clim_Lag", L), value = m[[paste0(v, "_z_25yr")]])
}
need <- paste0(rep(paste0(SHV, "_clim"), each = 29), "_Lag", 0:28)
ok <- complete.cases(D[, ..need])
cat(sprintf("climatology available for %.1f%% of rows | cases %d -> %d\n",
            100 * mean(ok), sum(D$outbreak_binary), sum(D$outbreak_binary[ok])))
D <- D[ok]
D[, `:=`(nc = sum(outbreak_binary), nr = sum(outbreak_binary == 0)), by = stratum_id]
D <- D[nc == 1L & nr > 0L][, c("nc", "nr") := NULL]
# z-score so the +0.5 SD contrast means the same thing as in the other models
for (v in paste0(SHV, "_clim")) for (L in 0:28) {
  cn <- paste0(v, "_Lag", L); set(D, j = cn, value = as.numeric(scale(D[[cn]])))
}
cat(sprintf("analytic sample after cleanup: %d cases, %.1f referents per case\n",
            sum(D$outbreak_binary), sum(D$outbreak_binary == 0) / sum(D$outbreak_binary)))

# 3. Three parameterisations on identical rows
SPECS <- list(
  "25-year climatological anomaly" = paste0(SHV, "_clim"),
  "90-day rolling anomaly"         = paste0(SHV, "_shock"),
  "Absolute conditions"            = SHV)

res <- rbindlist(lapply(names(SPECS), function(nm) {
  mv <- c(SPECS[[nm]], "anseriformes")
  av <- c(setNames(rep(list(NS2), length(SPECS[[nm]])), SPECS[[nm]]),
          list(anseriformes = LIN))
  o <- assemble_cc_dlnm(D, met_vars = mv, weekly_vars = character(0), max_lag = 28,
                        argvar_by_var = av, arglag = list(fun = "ns", df = 2))
  f <- suppressWarnings(fit_cc_inla_dlnm(o, fixed_prec = PREC))
  cat(sprintf("  %-32s WAIC %.1f\n", nm, f$waic$waic))
  base <- sub("_clim$|_shock$", "", SPECS[[nm]])
  rbindlist(lapply(seq_along(mv), function(i) {
    v <- mv[i]
    key <- if (v == "anseriformes") "anseriformes" else base[i]
    data.table(Parameterisation = nm, WAIC = round(f$waic$waic, 1),
               Predictor = unname(LABS[key]),
               `Lag 0-7 days`   = cell(met_lag_effect(f, o, v, 0:7)),
               `Lag 8-14 days`  = cell(met_lag_effect(f, o, v, 8:14)),
               `Lag 15-21 days` = cell(met_lag_effect(f, o, v, 15:21)),
               `Lag 22-28 days` = cell(met_lag_effect(f, o, v, 22:28)))
  }))
}))

saveRDS(res, paste0(objects_folder, "si_climatology_anomaly.RDS"))
cat("\n=== Minnesota: three exposure parameterisations, identical rows ===\n")
print(res)
