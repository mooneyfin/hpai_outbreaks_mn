############ LEAVE-ONE-STRATUM-OUT CROSS-VALIDATION, THE TWO PRIMARY MODELS ############
# Replaces 5-fold. With 167 strata, 5-fold trains on 80% of the data and so understates what the
# reported model can do; LOO trains on 166/167, which is essentially the model we report.
#
# The unit is the STRATUM, not the row. INLA's built-in CPO leaves out one observation at a
# time, which is not a hold-out here: the stratum intercept stays estimable from that stratum's
# remaining days, so the left-out day is still partly fitted. Dropping the whole stratum and
# refitting is the only honest version.
#
# Score: the log probability the model gives the day the outbreak actually happened, out of that
# stratum's candidate days. Chance is log(1/k) for a stratum with k days. The stratum intercept
# cancels in the conditional likelihood, so it never has to be estimated for the held-out
# stratum - which is what makes this possible at all.

rm(list = ls())
project.folder = paste0(print(here::here()), '/')
source(paste0(project.folder, 'create_folder_structure.R'))
source(paste0(functions.folder, 'script_initiate.R'))
suppressMessages({library(INLA); library(Matrix)})
INLA::inla.setOption(num.threads = 2)
source(paste0(functions.folder, 'inla_dlnm_helpers.R'))

PREC <- 4; NS2 <- list(fun = "ns", df = 2); LIN <- list(fun = "lin")
# exposures come from INLA_PRIMARY_* in inla_dlnm_helpers.R - never redefined locally
SHOCK <- INLA_PRIMARY_SHOCK
LEVEL <- INLA_PRIMARY_MET
MODELS <- list(`Absolute conditions` = LEVEL, `90-day shock` = SHOCK)
PANELS <- c("Minnesota", "Other flyway states", "Pooled")

D <- setDT(readRDS(paste0(objects_folder, "case_crossover_df_shock_post7.RDS")))
if (anyDuplicated(names(D))) D <- D[, which(!duplicated(names(D))), with = FALSE]
pick <- function(pn) {
  if (pn == "Minnesota") D[state == "Minnesota"]
  else if (pn == "Pooled") copy(D) else D[state != "Minnesota"]
}
assemble <- function(d, vars) {
  mv <- c(vars, "anseriformes")
  av <- c(setNames(rep(list(NS2), length(vars)), vars), list(anseriformes = LIN))
  assemble_cc_dlnm(d, met_vars = mv, weekly_vars = character(0), max_lag = 28,
                   argvar_by_var = av, arglag = list(fun = "ns", df = 2))
}
eta_fixed <- function(fit, obj) {
  cbn <- unlist(obj$cb_cols)
  b <- fit$summary.fixed[match(cbn, rownames(fit$summary.fixed)), "mean"]
  as.numeric(as.matrix(obj$data[, cbn, drop = FALSE]) %*% b)
}
score_one <- function(eta, y) eta[y == 1][1] - (max(eta) + log(sum(exp(eta - max(eta)))))

res <- rbindlist(lapply(PANELS, function(pn) {
  dat <- pick(pn); st <- unique(dat$stratum_id)
  cat(sprintf("%s: %d strata\n", pn, length(st)))
  rbindlist(lapply(names(MODELS), function(mn) {
    vars <- MODELS[[mn]]
    ll <- vapply(seq_along(st), function(i) {
      tr <- dat[stratum_id != st[i]]; te <- dat[stratum_id == st[i]]
      o_tr <- assemble(tr, vars); o_te <- assemble(te, vars)
      f_tr <- suppressWarnings(fit_cc_inla_dlnm(o_tr, fixed_prec = PREC))
      out <- try(score_one(eta_fixed(f_tr, o_te), o_te$data$outbreak_binary), silent = TRUE)
      if (inherits(out, "try-error")) NA_real_ else out
    }, numeric(1))
    cat(sprintf("  %-20s LOO score %.1f (%d/%d usable)\n", mn, sum(ll, na.rm = TRUE),
                sum(!is.na(ll)), length(ll)))
    data.table(Panel = pn, Model = mn, stratum = st, ll = ll)
  }))
}))
saveRDS(res, paste0(objects_folder, "si_loo_raw.RDS"))

chance <- D[, .(ll0 = -log(.N)), by = stratum_id]
summ <- rbindlist(lapply(PANELS, function(pn) {
  w <- dcast(res[Panel == pn], stratum ~ Model, value.var = "ll")
  w <- merge(w, chance, by.x = "stratum", by.y = "stratum_id")
  d <- w[["90-day shock"]] - w[["Absolute conditions"]]
  z <- mean(d, na.rm = TRUE) / (sd(d, na.rm = TRUE) / sqrt(sum(!is.na(d))))
  rbindlist(lapply(names(MODELS), function(mn) data.table(
    Panel = pn, Model = mn, Strata = as.character(nrow(w)),
    Score = sprintf("%.1f", sum(w[[mn]], na.rm = TRUE)),
    `vs chance` = sprintf("%+.1f", sum(w[[mn]], na.rm = TRUE) - sum(w$ll0, na.rm = TRUE)),
    `Anomaly - absolute, per stratum` = if (mn == "90-day shock")
      sprintf("%+.4f", mean(d, na.rm = TRUE)) else "reference",
    z = if (mn == "90-day shock") sprintf("%+.2f", z) else "-")))
}))
saveRDS(summ, paste0(objects_folder, "si_cv_primary.RDS"))
cat("\n=== leave-one-stratum-out, the two primary models ===\n"); print(summ)
