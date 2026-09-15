############ CLIMATOLOGICAL ANOMALY, ABSOLUTE UNITS — THE FIXED VERSION ############
# The standardised climatological anomaly z = (x - mean)/sd is broken for the skewed exposures.
# The sd is near zero in cell-weeks where runoff or precipitation is historically almost always
# zero, so dividing by it produces enormous values: 53,031 runoff and 67,559 precipitation
# cell-days with |z| > 10, against ZERO for temperature and soil moisture. The resulting
# exposure is dominated by a handful of absurd observations, and the null we previously reported
# for precipitation was fitted on that.
#
# Fix: use the ABSOLUTE anomaly, x - (25-year weekly mean), in natural units. No division by an
# unstable sd, well behaved for every exposure, and directly interpretable - "2 degC above the
# 25-year normal for this week here", "5 mm above normal".
#
# Weekly lag bins, because the climatology is weekly: every day inside a week is compared against
# the same normal, so a daily lag axis estimates structure the exposure does not have. Basis
# reselected by WAIC rather than carried over from the daily model.

rm(list = ls())
project.folder = paste0(print(here::here()), '/')
source(paste0(project.folder, 'create_folder_structure.R'))
source(paste0(functions.folder, 'script_initiate.R'))
suppressMessages({library(INLA); library(Matrix)})
INLA::inla.setOption(num.threads = 2)
source(paste0(functions.folder, 'inla_dlnm_helpers.R'))

PREC <- 4; NWEEK <- 4
VARS <- c("temperature", "precipitation", "soil_moisture", "wind_speed")
LAB <- c(setNames(c("Temperature", "Precipitation", "Soil moisture", "Wind speed"), VARS),
         anseriformes = "Anseriformes")
cell <- function(e) sprintf("%.2f (%.2f, %.2f)", e[["OR"]], e[["low"]], e[["high"]])

clim <- rbindlist(lapply(
  list.files(env_data, pattern = "flyway8_climatology_1997_2021_w", full.names = TRUE),
  fread, colClasses = list(character = "zone_id")))
daily <- setDT(readRDS(paste0(objects_folder, multistate_ts_daily_rds)))
daily[, date := as.Date(date)][, week_idx := pmin(as.integer(strftime(date, "%j")) %/% 7L, 51L)]
d <- merge(daily[, c("zone_id", "date", "week_idx", VARS), with = FALSE],
           clim[, c("zone_id", "week_idx", paste0(VARS, "_mean")), with = FALSE],
           by = c("zone_id", "week_idx"), all.x = TRUE)
for (v in VARS) set(d, j = paste0(v, "_cabs"), value = d[[v]] - d[[paste0(v, "_mean")]])
ABS <- paste0(VARS, "_cabs")
saveRDS(d[, c("zone_id", "date", ABS), with = FALSE],
        paste0(objects_folder, "daily_climatology_absolute.RDS"))

# attach 0-27 day lags of the absolute anomaly onto the case-crossover panel
cc <- setDT(readRDS(paste0(objects_folder, "case_crossover_df_climanom_post7.RDS")))
if (anyDuplicated(names(cc))) cc <- cc[, which(!duplicated(names(cc))), with = FALSE]
cc[, date := as.Date(date)]
src <- d[, c("zone_id", "date", ABS), with = FALSE]
for (L in 0:27) {
  s <- copy(src)[, date := date + L]
  setnames(s, ABS, paste0(ABS, "_Lag", L))
  cc <- merge(cc, s, by = c("zone_id", "date"), all.x = TRUE, sort = FALSE)
}
need <- paste0(rep(ABS, each = 28), "_Lag", 0:27)
cc <- cc[complete.cases(cc[, ..need])]
cc[, `:=`(k = sum(outbreak_binary), n = .N), by = stratum_id]
cc <- cc[k == 1L & n > 1L]
cat(sprintf("panel: %d cases, %d strata\n", sum(cc$outbreak_binary), uniqueN(cc$stratum_id)))

week_mat <- function(dat, v) sapply(seq_len(NWEEK) - 1L, function(w) {
  cols <- intersect(paste0(v, "_Lag", (w * 7):(w * 7 + 6)), names(dat))
  as.numeric(scale(rowMeans(dat[, ..cols], na.rm = TRUE)))
})
fit_w <- function(dat, av, al) {
  bases <- list(); dd <- data.frame(outbreak_binary = dat$outbreak_binary)
  for (v in c(ABS, "anseriformes")) {
    a <- if (v == "anseriformes") list(fun = "lin") else av
    cb <- dlnm::crossbasis(week_mat(dat, v), lag = c(0, NWEEK - 1), argvar = a, arglag = al)
    colnames(cb) <- paste0("cb_", v, ".", colnames(cb))
    bases[[v]] <- cb; dd <- cbind(dd, as.data.frame(unclass(cb)))
  }
  dd$id_stratum <- as.integer(as.factor(dat$stratum_id))
  obj <- list(data = dd, bases = bases, cb_cols = lapply(bases, colnames),
              met_vars = names(bases), contrast_steps = list(), weekly = list(),
              weekly_vars = character(0), max_lag = NWEEK - 1, max_lag_by_var = list(),
              argvar = av, argvar_by_var = list(), arglag = al, arglag_by_var = list())
  list(obj = obj, fit = suppressWarnings(fit_cc_inla_dlnm(obj, fixed_prec = PREC)))
}
shape <- function(s) if (s == "lin") list(fun = "lin") else list(fun = "ns", df = 2)
sel <- rbindlist(lapply(c("lin", "ns2"), function(vs) rbindlist(lapply(c("lin", "ns2"), function(ls) {
  r <- try(fit_w(cc, shape(vs), shape(ls)), silent = TRUE)
  if (inherits(r, "try-error")) return(NULL)
  cat(sprintf("  argvar %-3s arglag %-3s WAIC %.1f\n", vs, ls, r$fit$waic$waic))
  data.table(var_shape = vs, lag_shape = ls, WAIC = r$fit$waic$waic)
}))))
setorder(sel, WAIC); best <- sel[1]
cat(sprintf("selected: argvar %s, arglag %s (WAIC %.1f)\n\n",
            best$var_shape, best$lag_shape, best$WAIC))

PANELS <- c("Minnesota", "Other flyway states", "Pooled")
out <- rbindlist(lapply(PANELS, function(pn) {
  dd <- if (pn == "Minnesota") cc[state == "Minnesota"]
        else if (pn == "Pooled") copy(cc) else cc[state != "Minnesota"]
  r <- fit_w(dd, shape(best$var_shape), shape(best$lag_shape))
  rbindlist(lapply(names(r$obj$bases), function(v)
    rbindlist(lapply(seq_len(NWEEK) - 1L, function(w) data.table(
      Panel = pn, n = sum(dd$outbreak_binary),
      Exposure = unname(LAB[sub("_cabs$", "", v)]),
      `Lag week` = paste0("Week ", w),
      est = cell(met_lag_effect(r$fit, r$obj, v, w)))))))
}))
saveRDS(list(selection = sel, estimates = out),
        paste0(objects_folder, "si_climatology_absolute.RDS"))
print(dcast(out, Exposure + Panel ~ `Lag week`, value.var = "est"))
