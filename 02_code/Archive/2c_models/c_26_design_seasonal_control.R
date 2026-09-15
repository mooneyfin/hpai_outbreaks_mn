############ SUPPLEMENT — REFERENT DESIGN AND SEASONAL CONTROL ############
# Why the primary uses a fixed calendar partition rather than symmetric bidirectional
# referents, argued from a measurement rather than from estimator theory.
#
# Panel A. How much season each design leaves inside a stratum. The 25-year weekly normal is
#   pure season - no weather in it at all - so the spread of that normal WITHIN a stratum is
#   exactly the seasonal contrast a design still permits. Note this is not what a trend-offset
#   diagnostic measures: a symmetric design balances referents either side of the case, so the
#   MEAN case-minus-referent difference cancels by construction while the SPREAD stays wide.
#
# Panel B. The same model fitted under both designs. Meteorological associations appear where
#   the residual seasonal spread is large and vanish where it is small; waterfowl holds under
#   both.

rm(list = ls())
project.folder = paste0(print(here::here()), '/')
source(paste0(project.folder, 'create_folder_structure.R'))
source(paste0(functions.folder, 'script_initiate.R'))
suppressMessages({library(INLA); library(Matrix); library(flextable); library(officer)})
INLA::inla.setOption(num.threads = 1)
source(paste0(functions.folder, 'inla_dlnm_helpers.R'))

PREC <- 4; NS2 <- list(fun = "ns", df = 2); LIN <- list(fun = "lin")
VARS <- c("temperature", "precipitation", "soil_moisture", "runoff", "wind_speed")
LB   <- c(temperature = "Temperature", precipitation = "Precipitation",
          soil_moisture = "Soil moisture", runoff = "Runoff",
          wind_speed = "Wind speed", anseriformes = "Anseriformes")
cell <- function(e) sprintf("%.2f (%.2f, %.2f)", e[["OR"]], e[["low"]], e[["high"]])

DES <- c(`Symmetric bidirectional, +/-14 to 28 days` = "case_crossover_df_sym_g14.RDS",
         `Time-stratified month, no washout`         = "case_crossover_df_timestrat_month.RDS",
         `Time-stratified month, post-only 7 days`   = "case_crossover_df_timestrat_month_post7.RDS")

read_design <- function(f) {
  d <- setDT(readRDS(paste0(objects_folder, f)))
  if (anyDuplicated(names(d))) d <- d[, which(!duplicated(names(d))), with = FALSE]
  d[, date := as.Date(date)][]
}

# ---- Panel A: residual seasonal contrast ----
clim <- rbindlist(lapply(
  list.files(env_data, pattern = "flyway8_climatology_1997_2021_w", full.names = TRUE),
  fread, select = c("zone_id", "week_idx", "temperature_mean", "soil_moisture_mean")))
tot_t <- sd(clim$temperature_mean, na.rm = TRUE)
tot_s <- sd(clim$soil_moisture_mean, na.rm = TRUE)

panelA <- rbindlist(lapply(names(DES), function(dn) {
  d <- read_design(DES[[dn]])[, .(zone_id, date, stratum_id)]
  d[, week_idx := pmin(as.integer(strftime(date, "%j")) %/% 7L, 51L)]
  m <- merge(d, clim, by = c("zone_id", "week_idx"), all.x = TRUE)
  s <- m[, .(span = diff(range(as.integer(date))),
             t_sd = sd(temperature_mean, na.rm = TRUE),
             s_sd = sd(soil_moisture_mean, na.rm = TRUE)), by = stratum_id]
  data.table(Design = dn,
             `Stratum span (days)` = sprintf("%.1f", mean(s$span, na.rm = TRUE)),
             `Residual seasonal temperature (degC)` = sprintf("%.2f", mean(s$t_sd, na.rm = TRUE)),
             `% of annual seasonal range` = sprintf("%.0f%%",
               100 * mean(s$t_sd, na.rm = TRUE) / tot_t),
             `Residual seasonal soil moisture` = sprintf("%.4f", mean(s$s_sd, na.rm = TRUE)),
             `% of annual range` = sprintf("%.0f%%",
               100 * mean(s$s_sd, na.rm = TRUE) / tot_s))
}))

# ---- Panel B: same model, each design ----
panelB <- rbindlist(lapply(names(DES)[c(1, 3)], function(dn) {
  D <- read_design(DES[[dn]])
  rbindlist(lapply(c("Minnesota", "Pooled"), function(pn) {
    d  <- if (pn == "Minnesota") D[state == "Minnesota"] else copy(D)
    mv <- c(VARS, "anseriformes")
    av <- c(setNames(rep(list(NS2), length(VARS)), VARS), list(anseriformes = LIN))
    o <- assemble_cc_dlnm(d, met_vars = mv, weekly_vars = character(0), max_lag = 28,
                          argvar_by_var = av, arglag = list(fun = "ns", df = 2))
    f <- suppressWarnings(fit_cc_inla_dlnm(o, fixed_prec = PREC))
    rbindlist(lapply(mv, function(v) data.table(
      Design = dn, Panel = pn, Predictor = unname(LB[v]),
      `Lag 0-7 days`   = cell(met_lag_effect(f, o, v, 0:7)),
      `Lag 8-14 days`  = cell(met_lag_effect(f, o, v, 8:14)),
      `Lag 22-28 days` = cell(met_lag_effect(f, o, v, 22:28)))))
  }))
}))

saveRDS(list(seasonal_control = panelA, design_comparison = panelB),
        paste0(objects_folder, "tableS2_design.RDS"))

is_sig <- function(x) {
  m <- regmatches(x, regexec("\\(([0-9.]+), ([0-9.]+)\\)", x))
  vapply(m, function(z) length(z) == 3 && (as.numeric(z[2]) > 1 | as.numeric(z[3]) < 1), TRUE)
}
mk <- function(t, bold_cols = character(0)) {
  ft <- flextable(t) |> fontsize(size = 10, part = "all") |> autofit() |>
    align(align = "center", part = "all") |> theme_booktabs()
  for (j in bold_cols) {
    hit <- which(is_sig(t[[j]])); if (length(hit)) ft <- bold(ft, i = hit, j = j)
  }
  ft
}
read_docx() |>
  body_add_flextable(mk(panelA)) |> body_add_par("") |>
  body_add_flextable(mk(panelB, grep("^Lag", names(panelB), value = TRUE))) |>
  print(target = file.path(tables_main_folder, "TableS2_design_seasonal_control.docx"))

cat("=== Panel A: season left inside a stratum ===\n"); print(panelA)
cat("\n=== Panel B: same model, both designs ===\n"); print(panelB)
cat("\nwrote TableS2_design_seasonal_control.docx\n")
