# c_10_triangulation_cox_ts.R
# Triangulation of the case-crossover finding with two season-controlled designs, each on
# BOTH Minnesota and the 7-state multistate panel:
#   (A) Andersen-Gill recurrent-event Cox (between-cell spatial contrast at fixed calendar
#       time; season controlled via monthly baseline strata; poultry density adjusted) — mirrors c_02.
#   (B) Season-adjusted quasi-Poisson time-series DLNM on case cells (within-cell temporal;
#       season controlled via an explicit ns(week) trend).
# Both use the weekly lag 0-4 wk crossbasis (ns df=2) and the 7 primary exposures, per +0.5 SD.
# If climate stays weak and anseriformes persists across designs & regions, the anomaly
# case-crossover conclusion (birds real, climate seasonal) is corroborated by orthogonal methods.

rm(list = ls())
project.folder = paste0(print(here::here()), '/')
source(paste0(project.folder, 'create_folder_structure.R'))
source(paste0(functions.folder, 'script_initiate.R'))
for (fn in c("year", "month", "isoweek", "isoyear", "wday", "week"))
  conflict_prefer(fn, "lubridate", quiet = TRUE)
library(splines)

predictors  <- c("runoff", "soil_moisture", "temperature", "precipitation",
                 "wind_speed", "snow_cover", "anseriformes")
pred_labels <- c(runoff = "Runoff", soil_moisture = "Soil Moisture", temperature = "Temperature",
                 precipitation = "Precipitation", wind_speed = "Wind Speed",
                 snow_cover = "Snow Cover", anseriformes = "Anseriformes")
MAXW <- 4

# log1p the skewed exposures + z-scale all; poultry density (birds_au) → z (MN only in multistate).
prep_z <- function(w) {
  setDT(w)
  for (v in c("runoff", "precipitation", "anseriformes")) w[[v]] <- log1p(w[[v]])
  w[, birds_au_z := as.numeric(scale(log1p(fcoalesce(as.numeric(birds_au), 0))))]
  for (v in predictors) w[[v]] <- as.numeric(scale(w[[v]]))
  w[, month := lubridate::month(as.Date(week_start))]
  w[]
}

# MN weekly panel: aggregate the MN daily modeling panel to zone-weeks (mirrors c_02).
prep_mn <- function() {
  d <- setDT(readRDS(paste0(objects_folder, "fishnet_timeseries_daily_modeling.RDS")))
  if (max(d$runoff, na.rm = TRUE) < 1) d[, `:=`(runoff = runoff * 1000, precipitation = precipitation * 1000)]
  d[, `:=`(week_id = paste0(isoyear(date), "-W", formatC(isoweek(date), width = 2, flag = "0")),
           dow = wday(date, week_start = 1))]
  w <- d[, .(week_start = min(date), year = isoyear(min(date)),
             outbreak_binary = as.integer(any(outbreak_binary == 1)), birds_au = first(birds_au),
             temperature = mean(temperature, na.rm = TRUE), soil_moisture = mean(soil_moisture, na.rm = TRUE),
             snow_cover = mean(snow_cover, na.rm = TRUE), wind_speed = mean(wind_speed, na.rm = TRUE),
             anseriformes = { mv <- .SD[dow == 1, anseriformes]
                              if (length(mv) && !is.na(mv[1])) mv[1] else first(.SD$anseriformes) },
             precipitation = sum(precipitation, na.rm = TRUE), runoff = sum(runoff, na.rm = TRUE)),
         by = .(zone_id, week_id)]
  prep_z(w[year == 2022])
}

# Multistate weekly panel (already zone-week); restrict to 2022.
prep_ms <- function() {
  w <- setDT(readRDS(paste0(objects_folder, "multistate_fishnet_timeseries_weekly.RDS")))
  if (max(w$runoff, na.rm = TRUE) < 1) w[, `:=`(runoff = runoff * 1000, precipitation = precipitation * 1000)]
  w[, week_start := as.Date(week_start)]
  prep_z(w[year(week_start) == 2022])
}

# Cumulative effect (per +0.5 SD over lag 0-4 wk) for each predictor from a fitted model + crossbases.
cum_tab <- function(model, cbs) {
  beta <- coef(model); V <- vcov(model)
  rbindlist(lapply(predictors, function(v) {
    nm <- paste0("cb_", v); ix <- grep(paste0("^", nm, "v"), names(beta))
    b <- beta[ix]; Vs <- V[ix, ix, drop = FALSE]
    lb <- do.call(dlnm::onebasis, c(list(x = 0:MAXW), attr(cbs[[nm]], "arglag")))
    ct <- 0.5 * colSums(lb); lo <- sum(ct * b); se <- sqrt(as.numeric(t(ct) %*% Vs %*% ct))
    data.table(label = pred_labels[[v]],
               cell = sprintf("%.2f (%.2f,%.2f)%s", exp(lo), exp(lo - 1.96 * se), exp(lo + 1.96 * se),
                              ifelse(exp(lo - 1.96 * se) > 1 | exp(lo + 1.96 * se) < 1, "*", "")))
  }))
}

make_cbs <- function(w) setNames(lapply(predictors, function(v)
  crossbasis(w[[v]], lag = c(0, MAXW), argvar = list(fun = "lin"),
             arglag = list(fun = "ns", df = 2), group = w$zone_id)), paste0("cb_", predictors))

# (A) Andersen-Gill Cox: case cells, calendar-week intervals, monthly strata, poultry adjusted.
fit_cox <- function(w) {
  setorder(w, zone_id, week_start)
  w[, `:=`(tstart = seq_len(.N) - 1L, tstop = seq_len(.N), event = as.integer(outbreak_binary),
           smonth = month), by = zone_id]
  ag <- w[zone_id %in% w[event == 1, unique(zone_id)]]
  cbs <- make_cbs(ag); list2env(cbs, environment())
  pou <- if (var(ag$birds_au_z, na.rm = TRUE) > 0) "+ birds_au_z" else ""
  f <- as.formula(paste("Surv(tstart,tstop,event) ~", paste(names(cbs), collapse = " + "),
                        pou, "+ strata(smonth)"))
  cum_tab(coxph(f, data = ag, ties = "efron", cluster = zone_id), cbs)
}

# (B) Season-adjusted quasi-Poisson time-series DLNM on case cells (ns(week) seasonal trend).
fit_ts <- function(w) {
  setorder(w, zone_id, week_start)
  ts <- w[zone_id %in% w[outbreak_binary == 1, unique(zone_id)]]
  ts[, week_num := as.integer(as.factor(week_start))]
  cbs <- make_cbs(ts); list2env(cbs, environment())
  f <- as.formula(paste("outbreak_binary ~", paste(names(cbs), collapse = " + "), "+ ns(week_num, df = 6)"))
  cum_tab(glm(f, data = ts, family = quasipoisson), cbs)
}

mn <- prep_mn(); ms <- prep_ms()
res <- Reduce(function(a, b) merge(a, b, by = "label", all = TRUE), list(
  setnames(fit_cox(copy(mn)), "cell", "CoxAG_MN"),
  setnames(fit_cox(copy(ms)), "cell", "CoxAG_MS"),
  setnames(fit_ts(copy(mn)),  "cell", "TimeSeries_MN"),
  setnames(fit_ts(copy(ms)),  "cell", "TimeSeries_MS")))
res[, label := factor(label, levels = pred_labels[predictors])]; setorder(res, label)

cat("\n================ CONVERGENCE: cumulative effect (per +0.5 SD, lag 0-4 wk); * = CI excludes 1 =========\n")
cat("(A) Cox Andersen-Gill = HR (season controlled via monthly strata; poultry adjusted)\n")
cat("(B) Season-adjusted Poisson time-series = rate ratio (season controlled via ns(week))\n\n")
print(res)
saveRDS(res, paste0(objects_folder, "multistate_triangulation_cox_ts.RDS"))
