############ INPUTS FOR S2 PANEL A AND S7 ############
# Two things the supplement needs that aren't cached anywhere.
#
# 1. The UNIDIRECTIONAL design row for S2's ladder. It was the original submission's design and
#    is where the reviewer's concern is most visible, so the ladder is incomplete without it.
#    The panel stores absolute exposures only, so the 90-day shock columns get attached here the
#    same way the other design panels do.
#
# 2. S7 — the season interaction in S1's layout: Model x Exposure x Season rather than
#    Model x Exposure x Panel, so it can be read line for line against S1a.
#    SPRING vs OUTSIDE SPRING, not four seasons: case counts are Spring 112, Fall 44, Winter 8,
#    Summer 3, and a crossbasis with 10 parameters per exposure cannot be fitted on 3 cases.

rm(list = ls())
project.folder = paste0(print(here::here()), '/')
source(paste0(project.folder, 'create_folder_structure.R'))
source(paste0(functions.folder, 'script_initiate.R'))
suppressMessages({library(INLA); library(Matrix)})
INLA::inla.setOption(num.threads = 1)
source(paste0(functions.folder, 'inla_dlnm_helpers.R'))

PREC <- 4; NS2 <- list(fun = "ns", df = 2); LIN <- list(fun = "lin")
# exposures come from INLA_PRIMARY_* in inla_dlnm_helpers.R - never redefined locally
SHOCK <- INLA_PRIMARY_SHOCK
LEVEL <- INLA_PRIMARY_MET
MODELS <- list(`Absolute conditions` = LEVEL, `90-day shock` = SHOCK)
LAB <- INLA_PRIMARY_LAB
WIN <- INLA_LAG_WINDOWS
cell <- function(e) sprintf("%.2f (%.2f, %.2f)", e[["OR"]], e[["low"]], e[["high"]])

fit_it <- function(d, vars) {
  mv <- c(vars, "anseriformes")
  av <- c(setNames(rep(list(NS2), length(vars)), vars), list(anseriformes = LIN))
  o <- assemble_cc_dlnm(d, met_vars = mv, weekly_vars = character(0), max_lag = 28,
                        argvar_by_var = av, arglag = list(fun = "ns", df = 2))
  list(obj = o, fit = suppressWarnings(fit_cc_inla_dlnm(o, fixed_prec = PREC)), mv = mv)
}
win_row <- function(r, v) as.list(setNames(
  vapply(names(WIN), function(k) cell(met_lag_effect(r$fit, r$obj, v, seq(WIN[[k]][1], WIN[[k]][2]))),
         character(1)), names(WIN)))

# ---- 1. unidirectional design row ----
ex <- setDT(readRDS(paste0(objects_folder, "daily_shock_predictors.RDS")))
ex[, date := as.Date(date)]
SHC <- setdiff(names(ex), c("zone_id", "date"))
attach_shocks <- function(cc) {
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
uni <- setDT(readRDS(paste0(objects_folder, "case_crossover_df_multistate_daily.RDS")))
if (anyDuplicated(names(uni))) uni <- uni[, which(!duplicated(names(uni))), with = FALSE]
uni <- attach_shocks(uni)
cat(sprintf("unidirectional panel after attaching shocks: %d cases, %d strata, %.1f refs/case\n",
            sum(uni$outbreak_binary), uniqueN(uni$stratum_id),
            sum(uni$outbreak_binary == 0) / sum(uni$outbreak_binary)))

uni_row <- rbindlist(lapply(names(MODELS), function(mn) {
  r <- fit_it(uni, MODELS[[mn]])
  rbindlist(lapply(r$mv, function(v)
    cbind(data.table(Model = mn, Design = "Unidirectional, -7 to -28 d",
                     Exposure = unname(LAB[sub("_shock$", "", v)])), as.data.table(win_row(r, v)))))
}))
saveRDS(uni_row, paste0(objects_folder, "si_unidirectional_row.RDS"))

# ---- 2. S7 season interaction, S1 layout ----
D <- setDT(readRDS(paste0(objects_folder, "case_crossover_df_shock_post7.RDS")))
if (anyDuplicated(names(D))) D <- D[, which(!duplicated(names(D))), with = FALSE]
D[, spring := season == "Spring"]
# a stratum belongs to whichever season its CASE day falls in, so strata aren't split
D[, spring_stratum := any(spring[outbreak_binary == 1]), by = stratum_id]
cat(sprintf("spring strata %d, outside %d\n",
            uniqueN(D[spring_stratum == TRUE]$stratum_id),
            uniqueN(D[spring_stratum == FALSE]$stratum_id)))

s7 <- rbindlist(lapply(names(MODELS), function(mn)
  rbindlist(lapply(c(Spring = TRUE, `Outside spring` = FALSE), function(sp) {
    d <- D[spring_stratum == sp]
    r <- fit_it(d, MODELS[[mn]])
    rbindlist(lapply(r$mv, function(v)
      cbind(data.table(Model = mn, Season = if (sp) "Spring" else "Outside spring",
                       Exposure = unname(LAB[sub("_shock$", "", v)]),
                       Cases = sum(d$outbreak_binary)), as.data.table(win_row(r, v)))))
  }), idcol = NULL)))
saveRDS(s7, paste0(objects_folder, "si_season_s1layout.RDS"))
cat("\n=== S7 (spring vs outside), 90-day anomaly, key terms ===\n")
print(s7[Model == "90-day shock" & Exposure %in% c("Soil moisture", "Anseriformes")])
cat("\nwrote si_unidirectional_row.RDS and si_season_s1layout.RDS\n")
