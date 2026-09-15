############ POST-EXCLUSION LENGTH, ALL EXPOSURES ############
# The cached sweep only kept temperature, soil moisture and waterfowl, which makes the table
# look like the exposures were picked. Refit keeping every term, one-month strata only (S2a
# already covers stratum width), so the "nothing changes from 0 to 14 days" claim can be read
# across the whole exposure set rather than a chosen three.

rm(list = ls())
project.folder = paste0(print(here::here()), '/')
source(paste0(project.folder, 'create_folder_structure.R'))
source(paste0(functions.folder, 'script_initiate.R'))
suppressMessages({library(INLA); library(Matrix)})
INLA::inla.setOption(num.threads = 2)
source(paste0(functions.folder, 'inla_dlnm_helpers.R'))

PREC <- 4; NS2 <- list(fun = "ns", df = 2); LIN <- list(fun = "lin")
SHOCK <- c("temperature_shock", "precipitation_shock", "soil_moisture_shock", "wind_speed_shock")
LEVEL <- c("temperature", "precipitation", "soil_moisture", "wind_speed")
MODELS <- list(`Absolute conditions` = LEVEL, `90-day anomaly` = SHOCK)
LAB <- c(temperature = "Temperature", precipitation = "Precipitation",
         soil_moisture = "Soil moisture", wind_speed = "Wind speed", anseriformes = "Anseriformes")
WINDOW_OF <- c(Temperature = "Lag 0-7 days", Precipitation = "Lag 0-7 days",
               `Soil moisture` = "Lag 0-7 days", `Wind speed` = "Lag 0-7 days",
               Anseriformes = "Lag 22-28 days")

a06 <- readLines(paste0(data.prep.folder, "a_07_build_referent_panels.R"))
a06 <- a06[1:(grep("^for \\(dw in", a06)[1] - 1L)]
a06 <- a06[!grepl("^\\s*(rm\\(list|project\\.folder\\s*=|source\\()", a06)]
a06env <- new.env(parent = globalenv()); eval(parse(text = paste(a06, collapse = "\n")), a06env)
build_timestrat <- get("build_timestrat", a06env); daily <- get("daily", a06env)

ex <- setDT(readRDS(paste0(objects_folder, "daily_shock_predictors.RDS")))
ex[, date := as.Date(date)]
SHC <- setdiff(names(ex), c("zone_id", "date"))
attach_shocks <- function(cc) {
  if (anyDuplicated(names(cc))) cc <- cc[, which(!duplicated(names(cc))), with = FALSE]
  cc[, date := as.Date(date)]
  for (L in 0:28) {
    s <- ex[, c(.(zone_id = zone_id, date = date + L), .SD), .SDcols = SHC]
    setnames(s, SHC, paste0(SHC, "_Lag", L))
    cc <- merge(cc, s, by = c("zone_id", "date"), all.x = TRUE, sort = FALSE)
  }
  for (v in SHC) for (L in 0:28) {
    cn <- paste0(v, "_Lag", L); set(cc, j = cn, value = as.numeric(scale(cc[[cn]])))
  }
  need <- paste0(rep(SHC, each = 29), "_Lag", 0:28)
  cc <- cc[complete.cases(cc[, ..need])]
  cc[, `:=`(k = sum(outbreak_binary), n = .N), by = stratum_id]
  cc[k == 1L & n > 1L]
}

PANELS <- c("Minnesota", "Other flyway states", "Pooled")
pick <- function(d, pn) {
  if (pn == "Minnesota") d[state == "Minnesota"]
  else if (pn == "Pooled") copy(d) else d[state != "Minnesota"]
}
cell <- function(e) sprintf("%.2f (%.2f, %.2f)", e[["OR"]], e[["low"]], e[["high"]])

out <- rbindlist(lapply(c(0L, 3L, 5L, 7L, 10L, 14L), function(wo) {
  P <- attach_shocks(build_timestrat(daily, block_months = 1L, washout = wo,
                                     washout_type = "post"))
  rbindlist(lapply(names(MODELS), function(mn) rbindlist(lapply(PANELS, function(pn) {
    d <- pick(P, pn); vars <- MODELS[[mn]]
    mv <- c(vars, "anseriformes")
    av <- c(setNames(rep(list(NS2), length(vars)), vars), list(anseriformes = LIN))
    o <- assemble_cc_dlnm(d, met_vars = mv, weekly_vars = character(0), max_lag = 28,
                          argvar_by_var = av, arglag = list(fun = "ns", df = 2))
    f <- suppressWarnings(fit_cc_inla_dlnm(o, fixed_prec = PREC))
    cat(sprintf("  post %2d d | %-20s | %-20s refs/case %.1f\n", wo, mn, pn,
                sum(d$outbreak_binary == 0) / sum(d$outbreak_binary)))
    rbindlist(lapply(mv, function(v) {
      lb <- unname(LAB[sub("_shock$", "", v)]); w <- INLA_LAG_WINDOWS[[WINDOW_OF[[lb]]]]
      data.table(`Post exclusion (days)` = wo, Model = mn, Panel = pn, Exposure = lb,
                 Window = WINDOW_OF[[lb]],
                 `Referents per case` = sprintf("%.1f", sum(d$outbreak_binary == 0) /
                                                  sum(d$outbreak_binary)),
                 est = cell(met_lag_effect(f, o, v, seq(w[1], w[2]))))
    }))
  }))))
}))
saveRDS(out, paste0(objects_folder, "si_post_exclusion_all.RDS"))
cat("\nwrote si_post_exclusion_all.RDS\n")
