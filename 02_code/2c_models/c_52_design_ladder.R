############ REFERENT DESIGN LADDER — THE FITS BEHIND TABLE S2 ############
# Moved out of d_16_supplement_tables.R, which was a plotting script fitting models. Every
# referent design in the ladder gets both primary parameterisations, all six exposures, pooled
# panel. No dWAIC column on purpose: each design is a DIFFERENT set of observations, so the WAICs
# do not sit on the same data and differencing them just rewards whichever design deletes the
# most rows.
#
# The unidirectional row (the original submission's design) is built by c_33 - it needs its own
# panel. The negative-control column comes from c_30. d_20 stitches the three together.

project.folder = paste0(print(here::here()), '/')
src <- readLines(paste0(project.folder, 'create_folder_structure.R'))
eval(parse(text = paste(src[!grepl("^\\s*rm\\(list", src)], collapse = "\n")))
source(paste0(functions.folder, 'script_initiate.R'))
suppressMessages({library(INLA); library(Matrix)})
INLA::inla.setOption(num.threads = 4)
source(paste0(functions.folder, 'inla_dlnm_helpers.R'))

NS2 <- list(fun = "ns", df = 2); LIN <- list(fun = "lin")
# exposures come from INLA_PRIMARY_* in inla_dlnm_helpers.R - never redefined locally
MODELS <- list(`Absolute conditions` = INLA_PRIMARY_MET, `90-day shock` = INLA_PRIMARY_SHOCK)
BIRD   <- INLA_BIRD_VAR
cell <- function(e) sprintf("%.2f (%.2f, %.2f)", e[["OR"]], e[["low"]], e[["high"]])

# the selected design already carries its shock columns; the alternatives get them attached
D <- setDT(readRDS(paste0(objects_folder, "case_crossover_df_shock_post7.RDS")))
if (anyDuplicated(names(D))) D <- D[, which(!duplicated(names(D))), with = FALSE]
ex <- setDT(readRDS(paste0(objects_folder, "daily_shock_predictors.RDS")))
ex[, date := as.Date(date)]
SHOCK <- setdiff(names(ex), c("zone_id", "date"))
attach_shocks <- function(f) {
  cc <- setDT(readRDS(paste0(objects_folder, f)))
  if (anyDuplicated(names(cc))) cc <- cc[, which(!duplicated(names(cc))), with = FALSE]
  cc[, date := as.Date(date)]
  for (L in 0:28) {
    sr <- ex[, c(.(zone_id = zone_id, date = date + L), .SD), .SDcols = SHOCK]
    setnames(sr, SHOCK, paste0(SHOCK, "_Lag", L))
    cc <- merge(cc, sr, by = c("zone_id", "date"), all.x = TRUE, sort = FALSE)
  }
  for (v in SHOCK) for (L in 0:28) {
    cn <- paste0(v, "_Lag", L); set(cc, j = cn, value = as.numeric(scale(cc[[cn]])))
  }
  need <- paste0(rep(SHOCK, each = 29), "_Lag", 0:28)
  cc <- cc[complete.cases(cc[, ..need])]
  cc[, `:=`(nc = sum(outbreak_binary), nr = sum(outbreak_binary == 0)), by = stratum_id]
  cc[nc == 1L & nr > 0L]
}
DESIGNS <- list(
  `Symmetric bidirectional, +/-14 to 28 d`  = "case_crossover_df_sym_g14.RDS",
  `Time-stratified, two-month strata`       = "case_crossover_df_timestrat_bimonth_nowo.RDS",
  `One-month strata, 7-day post exclusion`  = NA)

# each exposure reported at its lag window: 0-7 d for weather, 22-28 d for waterfowl
fit_design <- function(d, vars, mn) {
  mv <- c(vars, BIRD)
  av <- c(setNames(rep(list(NS2), length(vars)), vars), setNames(list(LIN), BIRD))
  o <- assemble_cc_dlnm(d, met_vars = mv, weekly_vars = character(0), max_lag = 28,
                        argvar_by_var = av, arglag = NS2)
  f <- suppressWarnings(fit_cc_inla_dlnm(o, fixed_prec = INLA_PREC_PRIMARY))
  vals <- lapply(mv, function(v) {
    w <- if (v == BIRD) 22:28 else 0:7
    cell(met_lag_effect(f, o, v, w))
  })
  lab <- paste(INLA_PRIMARY_LAB[mv], ifelse(mv == BIRD, "22-28 d", "0-7 d"))
  as.list(setNames(vals, lab))
}

out <- rbindlist(lapply(names(DESIGNS), function(dn) {
  d <- if (is.na(DESIGNS[[dn]])) D else attach_shocks(DESIGNS[[dn]])
  refs <- sprintf("%.1f", sum(d$outbreak_binary == 0) / sum(d$outbreak_binary))
  cat(sprintf("%-42s cases %3d | %s refs/case\n", dn, sum(d$outbreak_binary), refs)); flush.console()
  rbindlist(lapply(names(MODELS), function(mn)
    cbind(data.table(Model = mn, Design = dn, Cases = sum(d$outbreak_binary),
                     `Referents per case` = refs),
          as.data.table(fit_design(d, MODELS[[mn]], mn)))))
}), fill = TRUE)
saveRDS(out, paste0(objects_folder, "si_design_ladder.RDS"))
cat("\nwrote si_design_ladder.RDS |", nrow(out), "rows\n"); print(out[, 1:5])
