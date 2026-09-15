############ CROSS-VALIDATION OF THE TWO PRIMARY MODELS ############
# The question: of the two parameterisations we report as primary - absolute conditions and the
# 90-day anomaly - which actually predicts better out of sample?
#
# How it works. Strata are split into 5 folds. Each fold is held out, the model is fitted on the
# other four, and we score the held-out strata by the log probability the model assigns to the
# day the outbreak actually happened. The stratum intercept cancels in the conditional
# likelihood, so it never has to be estimated for a stratum the model has not seen - which is
# what makes a case-crossover cross-validatable at all.
#
# Higher (less negative) is better. A model that knew nothing would spread probability evenly
# over a stratum's k days and score log(1/k); "vs chance" is the total gain on that.
#
# Comparisons are PAIRED on the stratum, because both models see the same fold assignment.
# Unpaired, fold noise swamps the difference at this sample size.

rm(list = ls())
project.folder = paste0(print(here::here()), '/')
source(paste0(project.folder, 'create_folder_structure.R'))
source(paste0(functions.folder, 'script_initiate.R'))
suppressMessages({library(INLA); library(Matrix)})
INLA::inla.setOption(num.threads = 2)
source(paste0(functions.folder, 'inla_dlnm_helpers.R'))

SEED <- 20260911; K <- 5; PREC <- 4
NS2 <- list(fun = "ns", df = 2); LIN <- list(fun = "lin")
SHOCK <- c("temperature_shock", "precipitation_shock", "soil_moisture_shock", "wind_speed_shock")
LEVEL <- c("temperature", "precipitation", "soil_moisture", "wind_speed")
MODELS <- list(`Absolute conditions` = LEVEL, `90-day anomaly` = SHOCK)
PANELS <- c("Minnesota", "Other flyway states", "Pooled")

D <- setDT(readRDS(paste0(objects_folder, "case_crossover_df_shock_post7.RDS")))
if (anyDuplicated(names(D))) D <- D[, which(!duplicated(names(D))), with = FALSE]
pick <- function(pn) {
  if (pn == "Minnesota") D[state == "Minnesota"]
  else if (pn == "Pooled") copy(D) else D[state != "Minnesota"]
}

# exposure terms only; the stratum intercept cancels in the conditional likelihood
eta_fixed <- function(fit, obj) {
  cbn <- unlist(obj$cb_cols)
  b <- fit$summary.fixed[match(cbn, rownames(fit$summary.fixed)), "mean"]
  as.numeric(as.matrix(obj$data[, cbn, drop = FALSE]) %*% b)
}
cond_ll <- function(eta, y, st) {
  d <- data.table(eta, y, st)
  d[, .(ll = eta[y == 1][1] - (max(eta) + log(sum(exp(eta - max(eta)))))), by = st]
}
assemble <- function(d, vars) {
  mv <- c(vars, "anseriformes")
  av <- c(setNames(rep(list(NS2), length(vars)), vars), list(anseriformes = LIN))
  assemble_cc_dlnm(d, met_vars = mv, weekly_vars = character(0), max_lag = 28,
                   argvar_by_var = av, arglag = list(fun = "ns", df = 2))
}

res <- rbindlist(lapply(PANELS, function(pn) {
  dat <- pick(pn)
  set.seed(SEED)
  st   <- unique(dat$stratum_id)
  fold <- setNames(sample(rep_len(1:K, length(st))), st)
  # chance baseline: one over the number of eligible days in each stratum
  chance <- dat[, .(ll0 = -log(.N)), by = stratum_id]

  per <- rbindlist(lapply(names(MODELS), function(mn) {
    ll <- rbindlist(lapply(1:K, function(k) {
      tr <- dat[stratum_id %in% names(fold)[fold != k]]
      te <- dat[stratum_id %in% names(fold)[fold == k]]
      o_tr <- assemble(tr, MODELS[[mn]]); o_te <- assemble(te, MODELS[[mn]])
      f_tr <- suppressWarnings(fit_cc_inla_dlnm(o_tr, fixed_prec = PREC))
      cond_ll(eta_fixed(f_tr, o_te), o_te$data$outbreak_binary, te$stratum_id)
    }))
    cat(sprintf("  %-20s %-20s score %.1f\n", pn, mn, sum(ll$ll)))
    data.table(Panel = pn, Model = mn, stratum = ll$st, ll = ll$ll)
  }))

  w <- dcast(per, Panel + stratum ~ Model, value.var = "ll")
  w <- merge(w, chance, by.x = "stratum", by.y = "stratum_id")
  d  <- w[["90-day anomaly"]] - w[["Absolute conditions"]]
  z  <- mean(d, na.rm = TRUE) / (sd(d, na.rm = TRUE) / sqrt(sum(!is.na(d))))
  rbindlist(lapply(names(MODELS), function(mn) data.table(
    Panel = pn, Model = mn,
    Strata = as.character(nrow(w)),
    Score = sprintf("%.1f", sum(w[[mn]], na.rm = TRUE)),
    `vs chance` = sprintf("%+.1f", sum(w[[mn]], na.rm = TRUE) - sum(w$ll0, na.rm = TRUE)),
    `Anomaly - absolute, per stratum` = if (mn == "90-day anomaly")
      sprintf("%+.4f", mean(d, na.rm = TRUE)) else "reference",
    z = if (mn == "90-day anomaly") sprintf("%+.2f", z) else "-")))
}))

saveRDS(res, paste0(objects_folder, "si_cv_primary.RDS"))
cat("\n=== cross-validation, the two primary models ===\n")
print(res)
