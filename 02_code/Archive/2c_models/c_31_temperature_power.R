############ CAN THE ONE-MONTH DESIGN SEE A TEMPERATURE EFFECT THE SIZE BIDIRECTIONAL REPORTS? ############
# Under symmetric bidirectional referents, absolute temperature at lag 0-7 reads 1.46 (1.14, 1.87).
# Under the selected one-month time-stratified design it is null. I had been explaining that as
# seasonal confounding in the bidirectional design - but the negative control (c_30) says
# bidirectional is NOT seasonally confounded, so that explanation is dead.
#
# Which leaves two possibilities: temperature really is null, or the one-month design is too
# narrow to detect it. A 28-day stratum contains far less temperature variation than a 55-day
# one, and the negative control's own interval was much wider under our design (0.90-1.42 vs
# 0.92-1.15), which is the same problem showing up directly.
#
# So: simulate outcomes under a truth that looks like what bidirectional found, and count how
# often the one-month design recovers it.
#
# This is design analysis in the Gelman & Carlin (2014) sense - power, type S (wrong sign) and
# type M (exaggeration) for a PRE-SPECIFIED external effect size. It is not post-hoc power
# computed from our own point estimate, which would be circular.
#
# Simulation is exact for this design: conditional on one case per stratum, the case day is
# multinomial with probabilities proportional to exp(x'beta). Fitted with clogit, which is exact
# here and fast enough for a thousand replicates.

rm(list = ls())
project.folder = paste0(print(here::here()), '/')
source(paste0(project.folder, 'create_folder_structure.R'))
source(paste0(functions.folder, 'script_initiate.R'))
suppressMessages({library(INLA); library(Matrix); library(survival)})
INLA::inla.setOption(num.threads = 1)
source(paste0(functions.folder, 'inla_dlnm_helpers.R'))

set.seed(20260911)
NREP   <- 1000
LEVEL  <- c("temperature", "precipitation", "soil_moisture", "wind_speed")
NS2    <- list(fun = "ns", df = 2)
MAXLAG <- 28

load_cc <- function(f) {
  d <- setDT(readRDS(paste0(objects_folder, f)))
  if (anyDuplicated(names(d))) d <- d[, which(!duplicated(names(d))), with = FALSE]
  d
}
assemble <- function(d) {
  av <- c(setNames(rep(list(NS2), length(LEVEL)), LEVEL), list(anseriformes = list(fun = "lin")))
  assemble_cc_dlnm(d, met_vars = c(LEVEL, "anseriformes"), weekly_vars = character(0),
                   max_lag = MAXLAG, argvar_by_var = av, arglag = list(fun = "ns", df = 2))
}

# 1. the truth we're testing against: temperature coefficients as bidirectional fitted them
o_bi <- assemble(load_cc("case_crossover_df_sym_g14.RDS"))
f_bi <- suppressWarnings(fit_cc_inla_dlnm(o_bi, fixed_prec = 4))
tcols <- o_bi$cb_cols[["temperature"]]
beta_shape <- f_bi$summary.fixed[tcols, "mean"]
e_bi <- met_lag_effect(f_bi, o_bi, "temperature", 0:7)
cat(sprintf("bidirectional temperature 0-7 d: %.2f (%.2f, %.2f)\n",
            e_bi[["OR"]], e_bi[["low"]], e_bi[["high"]]))

# 2. the design under test
o <- assemble(load_cc("case_crossover_df_timestrat_month_post7.RDS"))
# crossbasis columns live in o$bases, not o$data
X <- do.call(cbind, lapply(names(o$bases), function(v) unclass(o$bases[[v]])))
colnames(X) <- unlist(o$cb_cols)
tidx <- match(o$cb_cols[["temperature"]], colnames(X))
strat <- o$data$id_stratum

# contrast vector for the cumulative 0-7 d effect at +0.5 SD, in this design's own basis
cb <- o$bases[["temperature"]]
ob <- function(x, a) do.call(dlnm::onebasis, c(list(x = x), a))
v  <- as.numeric(ob(0.5, attr(cb, "argvar")) - ob(0, attr(cb, "argvar")))
W  <- colSums(ob(0:7, attr(cb, "arglag")), na.rm = TRUE)
cvec <- as.vector(t(outer(v, W)))

# scale the bidirectional shape so its 0-7 contrast equals the target RR in THIS basis
scale_to <- function(target_rr) {
  b <- beta_shape
  s <- log(target_rr) / sum(cvec * b)
  b * s
}

sim_once <- function(bt) {
  eta <- as.numeric(X[, tidx, drop = FALSE] %*% bt)   # only temperature has an effect
  y <- integer(length(eta))
  for (s in split(seq_along(eta), strat)) {
    p <- exp(eta[s] - max(eta[s])); y[s[sample.int(length(s), 1L, prob = p)]] <- 1L
  }
  dd <- data.table(y = y, strat = strat); dd <- cbind(dd, as.data.table(X))
  fm <- as.formula(paste("y ~", paste(colnames(X), collapse = " + "), "+ strata(strat)"))
  fit <- try(clogit(fm, data = dd), silent = TRUE)
  if (inherits(fit, "try-error")) return(c(NA, NA))
  b <- coef(fit)[colnames(X)[tidx]]
  V <- vcov(fit)[colnames(X)[tidx], colnames(X)[tidx], drop = FALSE]
  est <- sum(cvec * b); se <- sqrt(as.numeric(t(cvec) %*% V %*% cvec))
  c(est, se)
}

TARGETS <- c(1.20, 1.30, 1.46, 1.60, 1.80, 2.00)
res <- rbindlist(lapply(TARGETS, function(tr) {
  bt <- scale_to(tr)
  m <- vapply(seq_len(NREP), function(i) sim_once(bt), numeric(2))
  est <- m[1, ]; se <- m[2, ]; ok <- !is.na(est)
  z <- est[ok] / se[ok]
  sig <- abs(z) > 1.96
  data.table(`True RR` = sprintf("%.2f", tr),
             Power = sprintf("%.0f%%", 100 * mean(sig)),
             `Type S (wrong sign | significant)` =
               if (any(sig)) sprintf("%.1f%%", 100 * mean(sign(est[ok][sig]) != sign(log(tr)))) else "-",
             `Type M (exaggeration | significant)` =
               if (any(sig)) sprintf("%.2fx", mean(abs(est[ok][sig])) / abs(log(tr))) else "-",
             `Median RR recovered` = sprintf("%.2f", exp(median(est[ok]))))
}))
cat("\n=== one-month time-stratified design, absolute temperature, lag 0-7 d ===\n")
cat(sprintf("%d simulations per row, alpha = 0.05, %d strata\n\n", NREP, uniqueN(strat)))
print(res)
saveRDS(res, paste0(objects_folder, "si_temperature_power.RDS"))
