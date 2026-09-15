############ THREE OUTSTANDING SENSITIVITIES ############
# 1. Clustering of strata within cells - cluster bootstrap resampling CELLS, not strata
# 2. Shock window - the (recent, baseline) pair is a free parameter and was never probed
# 3. Season - two thirds of cases fall in the spring wave, so does anything hold within it
#
# All on the locked spec: 4 anomalies ns(2) + waterfowl linear, ns(2) lag, 0-28 d,
# N(0, 0.5^2), one-month strata with a 7-day post exclusion.

rm(list = ls())
project.folder = paste0(print(here::here()), '/')
source(paste0(project.folder, 'create_folder_structure.R'))
source(paste0(functions.folder, 'script_initiate.R'))
suppressMessages({library(INLA); library(Matrix)})
INLA::inla.setOption(num.threads = 1)
source(paste0(functions.folder, 'inla_dlnm_helpers.R'))

PREC <- 4; NS2 <- list(fun = "ns", df = 2); LIN <- list(fun = "lin")
SHV  <- c("temperature", "soil_moisture", "runoff", "wind_speed")
SH   <- paste0(SHV, "_shock")
MV   <- c(SH, "anseriformes")
AV   <- c(setNames(rep(list(NS2), length(SH)), SH), list(anseriformes = LIN))
LAB  <- c(setNames(c("Temperature", "Soil moisture", "Runoff", "Wind speed"), SH),
          anseriformes = "Anseriformes")
cell <- function(e) sprintf("%.2f (%.2f, %.2f)", e[["OR"]], e[["low"]], e[["high"]])

fit_it <- function(d, vars = MV, av = AV) {
  o <- assemble_cc_dlnm(d, met_vars = vars, weekly_vars = character(0), max_lag = 28,
                        argvar_by_var = av, arglag = list(fun = "ns", df = 2))
  list(obj = o, fit = suppressWarnings(fit_cc_inla_dlnm(o, fixed_prec = PREC)))
}
key2 <- function(r) c(
  soilM = met_lag_effect(r$fit, r$obj, "soil_moisture_shock", 0:7)[["OR"]],
  bird  = met_lag_effect(r$fit, r$obj, "anseriformes", 22:28)[["OR"]])

D <- setDT(readRDS(paste0(objects_folder, "case_crossover_df_shock_post7.RDS")))
if (anyDuplicated(names(D))) D <- D[, which(!duplicated(names(D))), with = FALSE]

############ 1. CLUSTER BOOTSTRAP BY CELL ############
# 27% of Minnesota cases sit in cells that contribute more than one stratum. The conditional
# likelihood treats those strata as independent; if a cell carries unmeasured shared risk it
# does not. Resample CELLS with replacement and see whether the spread of the bootstrap
# distribution is wider than the model-based interval.
NBOOT <- 200
boot_one <- function(d, seed) {
  set.seed(seed)
  cells <- unique(d$zone_id)
  pick  <- sample(cells, length(cells), replace = TRUE)
  # a cell drawn twice must contribute two DISTINCT strata, or they collapse into one
  bs <- rbindlist(lapply(seq_along(pick), function(i) {
    z <- d[zone_id == pick[i]]
    z <- copy(z)[, stratum_id := paste0(stratum_id, "_b", i)]
    z
  }))
  tryCatch(key2(fit_it(bs)), error = function(e) c(soilM = NA, bird = NA))
}
cat("cluster bootstrap by cell (", NBOOT, "draws )\n")
bs_res <- rbindlist(lapply(c("Minnesota", "Pooled"), function(pn) {
  d <- if (pn == "Minnesota") D[state == "Minnesota"] else copy(D)
  r0 <- fit_it(d)
  e_s <- met_lag_effect(r0$fit, r0$obj, "soil_moisture_shock", 0:7)
  e_b <- met_lag_effect(r0$fit, r0$obj, "anseriformes", 22:28)
  bo <- do.call(rbind, lapply(1:NBOOT, function(i) boot_one(d, 1000 + i)))
  bo <- bo[complete.cases(bo), , drop = FALSE]
  cat(sprintf("  %s: %d/%d draws converged\n", pn, nrow(bo), NBOOT))
  data.table(
    panel = pn,
    term = c("Soil moisture 0-7 days", "Anseriformes 22-28 days"),
    model_based = c(cell(e_s), cell(e_b)),
    bootstrap = c(sprintf("%.2f (%.2f, %.2f)", e_s[["OR"]],
                          quantile(bo[, "soilM"], .025), quantile(bo[, "soilM"], .975)),
                  sprintf("%.2f (%.2f, %.2f)", e_b[["OR"]],
                          quantile(bo[, "bird"], .025), quantile(bo[, "bird"], .975))),
    width_ratio = round(c(
      diff(quantile(bo[, "soilM"], c(.025, .975))) / (e_s[["high"]] - e_s[["low"]]),
      diff(quantile(bo[, "bird"],  c(.025, .975))) / (e_b[["high"]] - e_b[["low"]])), 2))
}))
cat("\n=== 1. Cluster bootstrap by cell ===\n"); print(bs_res)

############ 3. SEASON RESTRICTION ############
# run before the shock-window rebuild because it reuses the existing panel
D[, mth := as.integer(strftime(as.Date(date), "%m"))]
D[, case_mth := mth[outbreak_binary == 1][1], by = stratum_id]
seas_res <- rbindlist(lapply(c("Minnesota", "Pooled"), function(pn) {
  d0 <- if (pn == "Minnesota") D[state == "Minnesota"] else copy(D)
  rbindlist(lapply(list(All = NULL, `Spring (Mar-May)` = 3:5,
                        `Outside spring` = c(1:2, 6:12)), function(mm) {
    d <- if (is.null(mm)) d0 else d0[case_mth %in% mm]
    n <- sum(d$outbreak_binary)
    if (n < 20) return(NULL)
    r <- fit_it(d)
    data.table(panel = pn, n = n,
               soilM = cell(met_lag_effect(r$fit, r$obj, "soil_moisture_shock", 0:7)),
               bird  = cell(met_lag_effect(r$fit, r$obj, "anseriformes", 22:28)))
  }), idcol = "subset")
}))
cat("\n=== 3. Season restriction ===\n"); print(seas_res)

############ 2. SHOCK WINDOW ############
daily <- setDT(readRDS(paste0(objects_folder, multistate_ts_daily_rds)))
setorder(daily, zone_id, date); daily[, date := as.Date(date)]
for (v in c("precipitation", "runoff")) set(daily, j = v, value = daily[[v]] * 1000)
LOGV <- c("runoff")

build_shock_panel <- function(recent, baseline) {
  dd <- copy(daily)
  for (v in SHV) {
    x <- if (v %in% LOGV) log1p(pmax(dd[[v]], 0)) else dd[[v]]
    set(dd, j = "..x", value = x)
    dd[, (paste0(v, "_shock")) := frollmean(..x, recent, align = "right") -
         shift(frollmean(..x, baseline, align = "right"), 1), by = zone_id]
  }
  ex <- dd[, c("zone_id", "date", SH), with = FALSE]
  cc <- setDT(readRDS(paste0(objects_folder, "case_crossover_df_timestrat_month_post7.RDS")))
  if (anyDuplicated(names(cc))) cc <- cc[, which(!duplicated(names(cc))), with = FALSE]
  cc[, date := as.Date(date)]
  for (L in 0:28) {
    src <- ex[, c(.(zone_id = zone_id, date = date + L), .SD), .SDcols = SH]
    setnames(src, SH, paste0(SH, "_Lag", L))
    cc <- merge(cc, src, by = c("zone_id", "date"), all.x = TRUE, sort = FALSE)
  }
  for (v in SH) for (L in 0:28) {
    cn <- paste0(v, "_Lag", L); set(cc, j = cn, value = as.numeric(scale(cc[[cn]])))
  }
  need <- paste0(rep(SH, each = 29), "_Lag", 0:28)
  cc <- cc[complete.cases(cc[, ..need])]
  cc[, `:=`(nc = sum(outbreak_binary), nr = sum(outbreak_binary == 0)), by = stratum_id]
  cc[nc == 1L & nr > 0L]
}

WINDOWS <- list(c(7, 90), c(14, 90), c(7, 60), c(14, 60), c(7, 30), c(21, 90))
win_res <- rbindlist(lapply(WINDOWS, function(w) {
  P <- build_shock_panel(w[1], w[2])
  rbindlist(lapply(c("Minnesota", "Pooled"), function(pn) {
    d <- if (pn == "Minnesota") P[state == "Minnesota"] else copy(P)
    r <- fit_it(d)
    data.table(window = sprintf("%d-day recent vs %d-day baseline", w[1], w[2]),
               panel = pn, n = sum(d$outbreak_binary),
               WAIC = round(r$fit$waic$waic, 1),
               soilM = cell(met_lag_effect(r$fit, r$obj, "soil_moisture_shock", 0:7)),
               bird  = cell(met_lag_effect(r$fit, r$obj, "anseriformes", 22:28)))
  }))
}))
cat("\n=== 2. Shock window ===\n")
for (pn in c("Minnesota", "Pooled")) print(win_res[panel == pn][order(WAIC)])

saveRDS(list(cluster = bs_res, season = seas_res, window = win_res),
        paste0(objects_folder, "si_three_sensitivities.RDS"))
cat("\nsaved si_three_sensitivities.RDS\n")
