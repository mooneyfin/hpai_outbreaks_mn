############ FOLLOW-UPS: CLUSTER-ROBUST, BASELINE LENGTH, 25-YEAR CLIMATOLOGY ############
# 1. Cluster-robust intervals for every term, resampling CELLS rather than strata
# 2. Baseline length refit on COMMON rows - the earlier comparison was confounded, because a
#    shorter baseline needs less warm-up and so keeps more cases (167 vs 174 vs 175), and more
#    observations raise WAIC mechanically
# 3. A true climatological anomaly: (2022 weekly value - 1997-2021 weekly mean) / 25-yr sd.
#    Unlike the 90-day baseline this varies WITHIN a calendar month, so the stratum intercept
#    does not simply absorb it. Weekly and Minnesota-only, which is all the climatology covers.

rm(list = ls())
project.folder = paste0(print(here::here()), '/')
source(paste0(project.folder, 'create_folder_structure.R'))
source(paste0(functions.folder, 'script_initiate.R'))
suppressMessages({library(INLA); library(Matrix)})
INLA::inla.setOption(num.threads = 1)
source(paste0(functions.folder, 'inla_dlnm_helpers.R'))

PREC <- 4; NS2 <- list(fun = "ns", df = 2); LIN <- list(fun = "lin")
# exposures come from INLA_PRIMARY_* in inla_dlnm_helpers.R - never redefined locally
SHV  <- INLA_PRIMARY_MET
SH   <- paste0(SHV, "_shock")
MV   <- c(SH, "anseriformes")
AV   <- c(setNames(rep(list(NS2), length(SH)), SH), list(anseriformes = LIN))
WIN  <- INLA_LAG_WINDOWS
cell <- function(e) sprintf("%.2f (%.2f, %.2f)", e[["OR"]], e[["low"]], e[["high"]])

fit_it <- function(d, vars = MV, av = AV, ml = 28) {
  o <- assemble_cc_dlnm(d, met_vars = vars, weekly_vars = character(0), max_lag = ml,
                        argvar_by_var = av, arglag = list(fun = "ns", df = 2))
  list(obj = o, fit = suppressWarnings(fit_cc_inla_dlnm(o, fixed_prec = PREC)))
}
all_or <- function(r, vars) {
  unlist(lapply(vars, function(v) setNames(
    vapply(names(WIN), function(k) met_lag_effect(r$fit, r$obj, v, seq(WIN[[k]][1], WIN[[k]][2]))[["OR"]],
           numeric(1)), paste(v, names(WIN), sep = "|"))))
}

D <- setDT(readRDS(paste0(objects_folder, "case_crossover_df_shock_post7.RDS")))
if (anyDuplicated(names(D))) D <- D[, which(!duplicated(names(D))), with = FALSE]

############ 1. CLUSTER-ROBUST INTERVALS, EVERY TERM ############
NBOOT <- 200
boot_draw <- function(d, seed) {
  set.seed(seed)
  cells <- unique(d$zone_id)
  pick  <- sample(cells, length(cells), replace = TRUE)
  bs <- rbindlist(lapply(seq_along(pick), function(i)
    copy(d[zone_id == pick[i]])[, stratum_id := paste0(stratum_id, "_b", i)]))
  tryCatch(all_or(fit_it(bs), MV), error = function(e) NULL)
}
cat("cluster bootstrap,", NBOOT, "draws per panel\n")
clus <- rbindlist(lapply(c("Minnesota", "Pooled"), function(pn) {
  d  <- if (pn == "Minnesota") D[state == "Minnesota"] else copy(D)
  r0 <- fit_it(d)
  bo <- do.call(rbind, Filter(Negate(is.null), lapply(1:NBOOT, function(i) boot_draw(d, 2000 + i))))
  cat(sprintf("  %s: %d draws\n", pn, nrow(bo)))
  rbindlist(lapply(MV, function(v) rbindlist(lapply(names(WIN), function(k) {
    e  <- met_lag_effect(r0$fit, r0$obj, v, seq(WIN[[k]][1], WIN[[k]][2]))
    bq <- quantile(bo[, paste(v, k, sep = "|")], c(.025, .975), na.rm = TRUE)
    data.table(panel = pn, Predictor = unname(INLA_LABELS[sub("_shock$", "", v)]),
               window = k, model_based = cell(e),
               cluster_robust = sprintf("%.2f (%.2f, %.2f)", e[["OR"]], bq[1], bq[2]),
               width_ratio = round(unname(diff(bq)) / (e[["high"]] - e[["low"]]), 2))
  }))))
}))
cat("\n=== 1. Cluster-robust intervals ===\n")
print(clus[Predictor %in% c("Soil moisture", "Anseriformes")])

############ 2. BASELINE LENGTH ON COMMON ROWS ############
daily <- setDT(readRDS(paste0(objects_folder, multistate_ts_daily_rds)))
setorder(daily, zone_id, date); daily[, date := as.Date(date)]
for (v in c("precipitation", "runoff")) set(daily, j = v, value = daily[[v]] * 1000)

panel_for <- function(recent, baseline) {
  dd <- copy(daily)
  for (v in SHV) {
    x <- if (v == "runoff") log1p(pmax(dd[[v]], 0)) else dd[[v]]
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
  cc
}
KEEP <- D[, paste(zone_id, date)]                       # the 90-day panel's rows
base_res <- rbindlist(lapply(c(90, 60, 30), function(bl) {
  P <- panel_for(7, bl)
  P <- P[paste(zone_id, date) %in% KEEP]                # common rows, so WAIC is comparable
  for (v in SH) for (L in 0:28) {
    cn <- paste0(v, "_Lag", L); set(P, j = cn, value = as.numeric(scale(P[[cn]])))
  }
  rbindlist(lapply(c("Minnesota", "Pooled"), function(pn) {
    d <- if (pn == "Minnesota") P[state == "Minnesota"] else copy(P)
    r <- fit_it(d)
    data.table(baseline = sprintf("%d-day", bl), panel = pn, n = sum(d$outbreak_binary),
               WAIC = round(r$fit$waic$waic, 1),
               soilM = cell(met_lag_effect(r$fit, r$obj, "soil_moisture_shock", 0:7)),
               bird = cell(met_lag_effect(r$fit, r$obj, "anseriformes", 22:28)))
  }))
}))
cat("\n=== 2. Baseline length, common rows ===\n"); print(base_res[order(panel, WAIC)])

############ 3. 25-YEAR CLIMATOLOGICAL ANOMALY (MN, weekly) ############
an <- setDT(readRDS(paste0(objects_folder, "climate_anomaly_panel.RDS")))
ZC <- paste0(SHV, "_z_25yr")
an <- an[, c("zone_id", "week_start", ZC), with = FALSE]
an[, week_start := as.Date(week_start)]

MN <- D[state == "Minnesota"]
MN[, wk := as.Date(cut(as.Date(date), "week", start.on.monday = TRUE))]
# the climatology weeks may start on a different weekday, so join on nearest preceding week
setkey(an, zone_id, week_start)
mnz <- an[.(MN$zone_id, MN$wk), roll = TRUE, on = .(zone_id, week_start)]
for (v in ZC) set(MN, j = paste0(v, "_v"), value = mnz[[v]])

cat(sprintf("\nclimatology matched for %.0f%% of rows\n",
            100 * mean(!is.na(MN$temperature_z_25yr_v))))
# within-stratum correlation with the raw level: if ~1 the anomaly adds nothing here
for (v in SHV) {
  a <- MN[[paste0(v, "_z_25yr_v")]]; b <- MN[[paste0(v, "_Lag0")]]
  ok <- is.finite(a) & is.finite(b)
  w <- data.table(s = MN$stratum_id[ok], a = a[ok], b = b[ok])
  w[, `:=`(a = a - mean(a), b = b - mean(b)), by = s]
  cat(sprintf("  %-14s within-stratum cor(25-yr anomaly, level) = %+.3f\n", v, cor(w$a, w$b)))
}

saveRDS(list(cluster = clus, baseline = base_res),
        paste0(objects_folder, "si_followup_sensitivities.RDS"))
cat("\nsaved si_followup_sensitivities.RDS\n")
