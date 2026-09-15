############ HEAD-TO-HEAD: THREE EXPOSURE PARAMETERISATIONS ############
#   absolute conditions | 90-day rolling anomaly | 25-year climatological anomaly
#
# Two comparisons, because the panels are not the same size. The 90-day anomaly needs a
# 90-day warm-up and the daily series only starts 2021-11-15, so it loses 8 cases; the
# climatological anomaly needs no warm-up at all and keeps all 175. Fitting on the common
# rows makes WAIC comparable; fitting the climatology on its own full panel shows what it
# buys. (Exporting ERA5 back to ~2021-08-17 would remove the gap entirely.)

rm(list = ls())
project.folder = paste0(print(here::here()), '/')
source(paste0(project.folder, 'create_folder_structure.R'))
source(paste0(functions.folder, 'script_initiate.R'))
suppressMessages({library(INLA); library(Matrix)})
INLA::inla.setOption(num.threads = 1)
source(paste0(functions.folder, 'inla_dlnm_helpers.R'))

PREC <- 4; NS2 <- list(fun = "ns", df = 2); LIN <- list(fun = "lin")
# exposures come from INLA_PRIMARY_* in inla_dlnm_helpers.R - never redefined locally
VARS <- INLA_PRIMARY_MET
LABEL <- INLA_PRIMARY_LAB
cell <- function(e) sprintf("%.2f (%.2f, %.2f)", e[["OR"]], e[["low"]], e[["high"]])

# the climatology panel is the wider one (175 cases); bring the 90-day shock columns onto it
climp <- setDT(readRDS(paste0(objects_folder, "case_crossover_df_climanom_post7.RDS")))
if (anyDuplicated(names(climp))) climp <- climp[, which(!duplicated(names(climp))), with = FALSE]
shockp <- setDT(readRDS(paste0(objects_folder, "case_crossover_df_shock_post7.RDS")))
if (anyDuplicated(names(shockp))) shockp <- shockp[, which(!duplicated(names(shockp))), with = FALSE]
SHCOLS <- grep("_shock_Lag", names(shockp), value = TRUE)
climp[, date := as.Date(date)]; shockp[, date := as.Date(date)]
dat <- merge(climp, shockp[, c("zone_id", "date", SHCOLS), with = FALSE],
             by = c("zone_id", "date"), all.x = TRUE, sort = FALSE)

common <- complete.cases(dat[, ..SHCOLS])
cat(sprintf("climatology panel %d cases | rows with 90-day shocks too: %d cases\n",
            sum(dat$outbreak_binary), sum(dat$outbreak_binary[common])))

SPECS <- list(`Absolute conditions` = VARS,
              `90-day rolling anomaly` = paste0(VARS, "_shock"),
              `25-year climatological anomaly` = paste0(VARS, "_clim"))

fit_one <- function(d, exposures) {
  mv <- c(exposures, "anseriformes")
  av <- c(setNames(rep(list(NS2), length(exposures)), exposures), list(anseriformes = LIN))
  o <- assemble_cc_dlnm(d, met_vars = mv, weekly_vars = character(0), max_lag = 28,
                        argvar_by_var = av, arglag = list(fun = "ns", df = 2))
  list(obj = o, fit = suppressWarnings(fit_cc_inla_dlnm(o, fixed_prec = PREC)), mv = mv)
}
tidy <- function(r, label, panel, n) {
  rbindlist(lapply(r$mv, function(v) {
    base <- sub("_shock$|_clim$", "", v)
    data.table(Parameterisation = label, panel = panel, n = n,
               WAIC = round(r$fit$waic$waic, 1), Predictor = unname(LABEL[base]),
               `Lag 0-7 days`   = cell(met_lag_effect(r$fit, r$obj, v, 0:7)),
               `Lag 8-14 days`  = cell(met_lag_effect(r$fit, r$obj, v, 8:14)),
               `Lag 15-21 days` = cell(met_lag_effect(r$fit, r$obj, v, 15:21)),
               `Lag 22-28 days` = cell(met_lag_effect(r$fit, r$obj, v, 22:28)))
  }))
}

PANELS <- list(Minnesota = function(d) d[state == "Minnesota"],
               `Other flyway states` = function(d) d[state != "Minnesota"],
               Pooled = function(d) copy(d))

# (a) common rows - WAIC comparable across all three
cat("\n--- common rows ---\n")
common_res <- rbindlist(lapply(names(PANELS), function(pn) {
  d <- PANELS[[pn]](dat[common])
  rbindlist(lapply(names(SPECS), function(lab) {
    r <- fit_one(d, SPECS[[lab]])
    cat(sprintf("  %-13s %-32s n=%3d WAIC %.1f\n", pn, lab, sum(d$outbreak_binary),
                r$fit$waic$waic))
    tidy(r, lab, pn, sum(d$outbreak_binary))
  }))
}))

# (b) the climatological anomaly on its own full panel
cat("\n--- climatology on its full panel ---\n")
full_res <- rbindlist(lapply(names(PANELS), function(pn) {
  d <- PANELS[[pn]](dat)
  r <- fit_one(d, SPECS[["25-year climatological anomaly"]])
  cat(sprintf("  %-13s n=%3d WAIC %.1f\n", pn, sum(d$outbreak_binary), r$fit$waic$waic))
  tidy(r, "25-year climatological anomaly (full panel)", pn, sum(d$outbreak_binary))
}))

out <- rbind(common_res, full_res)
saveRDS(out, paste0(objects_folder, "si_three_parameterisations.RDS"))
fwrite(out, file.path(tables_main_folder, "TableS8_three_parameterisations.csv"))

cat("\n=== WAIC, common rows ===\n")
print(dcast(common_res[Predictor == "Anseriformes"], panel + n ~ Parameterisation,
            value.var = "WAIC"))
cat("\n=== key terms ===\n")
print(out[Predictor %in% c("Soil moisture", "Temperature", "Anseriformes"),
          .(Parameterisation, panel, Predictor, `Lag 0-7 days`, `Lag 22-28 days`)])
