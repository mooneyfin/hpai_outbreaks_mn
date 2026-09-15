############ CLUSTER-ROBUST INTERVALS ON THE CURRENT SPEC ############
# The cached cluster bootstrap (c_21 / si_followup_sensitivities.RDS) was run on the superseded
# exposure set - it carries runoff and no precipitation - so it can't be quoted next to the
# current estimates. Rerun here for both primary models on the pooled panel, which is where S1a
# wants the column.
#
# Resample CELLS, not cell-days. Strata within a cell are not independent: the same premises,
# the same land cover, the same local water, and repeat detections cluster inside 14 days. The
# model-based interval assumes independent strata and is consequently too narrow.

rm(list = ls())
project.folder = paste0(print(here::here()), '/')
source(paste0(project.folder, 'create_folder_structure.R'))
source(paste0(functions.folder, 'script_initiate.R'))
suppressMessages({library(INLA); library(Matrix)})
INLA::inla.setOption(num.threads = 1)
source(paste0(functions.folder, 'inla_dlnm_helpers.R'))

NBOOT <- 200
PREC  <- 4; NS2 <- list(fun = "ns", df = 2); LIN <- list(fun = "lin")
SHOCK <- c("temperature_shock", "precipitation_shock", "soil_moisture_shock", "wind_speed_shock")
LEVEL <- c("temperature", "precipitation", "soil_moisture", "wind_speed")
MODELS <- list(`Absolute` = LEVEL, `90-day anomaly` = SHOCK)
LAB <- c(temperature = "Temperature", precipitation = "Precipitation",
         soil_moisture = "Soil moisture", wind_speed = "Wind speed",
         anseriformes = "Anseriformes")
WIN <- INLA_LAG_WINDOWS

D <- setDT(readRDS(paste0(objects_folder, "case_crossover_df_shock_post7.RDS")))
if (anyDuplicated(names(D))) D <- D[, which(!duplicated(names(D))), with = FALSE]

fit_it <- function(d, vars) {
  mv <- c(vars, "anseriformes")
  av <- c(setNames(rep(list(NS2), length(vars)), vars), list(anseriformes = LIN))
  o <- assemble_cc_dlnm(d, met_vars = mv, weekly_vars = character(0), max_lag = 28,
                        argvar_by_var = av, arglag = list(fun = "ns", df = 2))
  list(obj = o, fit = suppressWarnings(fit_cc_inla_dlnm(o, fixed_prec = PREC)), mv = mv)
}
all_or <- function(r) unlist(lapply(r$mv, function(v)
  setNames(vapply(names(WIN), function(k)
    met_lag_effect(r$fit, r$obj, v, seq(WIN[[k]][1], WIN[[k]][2]))[["OR"]], numeric(1)),
    paste(v, names(WIN), sep = "|"))))

PANELS <- c("Minnesota", "Other flyway states", "Pooled")
pick_panel <- function(pn) {
  if (pn == "Minnesota") D[state == "Minnesota"]
  else if (pn == "Pooled") copy(D) else D[state != "Minnesota"]
}

# all three panels now: the cluster-robust interval is the reported one, so every row needs it
out <- rbindlist(lapply(PANELS, function(pn) rbindlist(lapply(names(MODELS), function(mn) {
  vars <- MODELS[[mn]]
  Dp <- pick_panel(pn)
  r0 <- fit_it(Dp, vars)
  cells <- unique(Dp$zone_id)
  bo <- do.call(rbind, Filter(Negate(is.null), lapply(1:NBOOT, function(i) {
    set.seed(3000 + i)
    pick <- sample(cells, length(cells), replace = TRUE)
    bs <- rbindlist(lapply(seq_along(pick), function(j)
      copy(Dp[zone_id == pick[j]])[, stratum_id := paste0(stratum_id, "_b", j)]))
    tryCatch(all_or(fit_it(bs, vars)), error = function(e) NULL)
  })))
  cat(sprintf("  %-20s %-16s %d usable draws\n", pn, mn, nrow(bo)))
  rbindlist(lapply(r0$mv, function(v) rbindlist(lapply(names(WIN), function(k) {
    e  <- met_lag_effect(r0$fit, r0$obj, v, seq(WIN[[k]][1], WIN[[k]][2]))
    bq <- quantile(bo[, paste(v, k, sep = "|")], c(.025, .975), na.rm = TRUE)
    data.table(Panel = pn, Model = mn, Predictor = unname(LAB[sub("_shock$", "", v)]), window = k,
               model_based = sprintf("%.2f (%.2f, %.2f)", e[["OR"]], e[["low"]], e[["high"]]),
               cluster_robust = sprintf("%.2f (%.2f, %.2f)", e[["OR"]], bq[1], bq[2]),
               width_ratio = round(unname(diff(bq)) / (e[["high"]] - e[["low"]]), 2))
  }))))
})))) 

saveRDS(out, paste0(objects_folder, "si_cluster_robust_current.RDS"))
cat("\n=== cluster-robust, all panels, current spec ===\n")
print(out[Predictor %in% c("Soil moisture", "Anseriformes")])
cat("\nwidth ratio > 1 means the model-based interval is too narrow\n")
