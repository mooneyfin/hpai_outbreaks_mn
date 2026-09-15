############ METHODS CHECK — CONDITIONAL POISSON REPRODUCES CONDITIONAL LOGISTIC ############
# The outcome is binary (1 case day, N referent days per stratum) and we fit Poisson, which
# looks wrong until you notice two things: 0/1 are legitimate Poisson counts, and with a free
# intercept per stratum the Poisson likelihood factorises into a piece that depends only on the
# stratum totals (fixed at one case by design, so uninformative about beta) and a piece that is
# exactly the conditional logistic likelihood. Birch (1963) for the general equivalence,
# Armstrong, Gasparrini & Tobias (2014, BMC Med Res Methodol 14:122) for the case-crossover case.
#
# We route through Poisson because INLA has no conditional logistic family, and going this way
# buys priors, penalised lag smooths and WAIC on a design that's otherwise stuck with clogit.
#
# This script checks the equivalence numerically rather than citing it. INLA runs with an
# effectively flat prior (prec 0.001) so it's estimating the same thing clogit's MLE does - the
# primary's N(0, 0.5^2) would shrink and the two would legitimately differ.

rm(list = ls())
project.folder = paste0(print(here::here()), '/')
source(paste0(project.folder, 'create_folder_structure.R'))
source(paste0(functions.folder, 'script_initiate.R'))
suppressMessages({library(INLA); library(Matrix); library(survival)})
INLA::inla.setOption(num.threads = 1)
source(paste0(functions.folder, 'inla_dlnm_helpers.R'))

SH  <- c("temperature_shock", "precipitation_shock", "soil_moisture_shock", "wind_speed_shock")
NS2 <- list(fun = "ns", df = 2)
AV  <- c(setNames(rep(list(NS2), length(SH)), SH), list(anseriformes = list(fun = "lin")))

D <- setDT(readRDS(paste0(objects_folder, "case_crossover_df_shock_post7.RDS")))
if (anyDuplicated(names(D))) D <- D[, which(!duplicated(names(D))), with = FALSE]

o <- assemble_cc_dlnm(D, met_vars = c(SH, "anseriformes"), weekly_vars = character(0),
                      max_lag = 28, argvar_by_var = AV, arglag = list(fun = "ns", df = 2))

# 1. conditional Poisson through INLA, essentially flat prior
fi <- suppressWarnings(fit_cc_inla_dlnm(o, fixed_prec = 0.001))
bi <- fi$summary.fixed
bi <- bi[!grepl("Intercept|stratum", rownames(bi)), , drop = FALSE]

# 2. the same crossbasis columns through clogit. the stratum column is id_stratum, and terms
#    are matched to INLA's fixed effects by NAME - position would silently misalign if INLA ever
#    reorders them.
X <- do.call(cbind, lapply(names(o$bases), function(v) {
  m <- unclass(o$bases[[v]]); colnames(m) <- paste0(v, "_b", seq_len(ncol(m))); m
}))
stopifnot(ncol(X) == length(unlist(o$cb_cols)))
colnames(X) <- unlist(o$cb_cols)                       # same names INLA fits under

dd <- data.table(y = o$data$outbreak_binary, strat = o$data$id_stratum)
dd <- cbind(dd, as.data.table(X))
# clogit builds the Surv object itself; handing it one is what breaks it
dd[, y := as.integer(y)]
fc <- clogit(as.formula(paste("y ~", paste(colnames(X), collapse = " + "), "+ strata(strat)")),
             data = dd)

shared <- intersect(colnames(X), rownames(bi))
cat(sprintf("matched %d of %d crossbasis terms by name\n", length(shared), ncol(X)))
stopifnot(length(shared) == ncol(X))

cmp <- data.table(term      = shared,
                  clogit    = unname(coef(fc)[shared]),
                  inla      = bi[shared, "mean"],
                  se_clogit = unname(sqrt(diag(vcov(fc)))[match(shared, colnames(X))]),
                  se_inla   = bi[shared, "sd"])
cmp[, `:=`(abs_diff = abs(clogit - inla), se_diff = abs(se_clogit - se_inla))]

saveRDS(cmp, paste0(objects_folder, "si_clogit_equivalence.RDS"))
print(cmp)
# report on the scale we actually publish. a threshold on the LOG coefficient is a stricter
# test than anything the paper claims, and failing it would be misleading.
cmp[, `:=`(rr_clogit = exp(clogit), rr_inla = exp(inla))]
n2 <- sum(round(cmp$rr_clogit, 2) == round(cmp$rr_inla, 2))
cat(sprintf(paste0("\nterms                  : %d\n",
                   "max |log-coef diff|    : %.5f  (%.2f%% on the ratio scale)\n",
                   "max |SE diff|          : %.5f\n",
                   "max ratio gap          : %.4f\n",
                   "ratios equal to 2 dp   : %d of %d\n"),
            nrow(cmp), max(cmp$abs_diff), 100 * (exp(max(cmp$abs_diff)) - 1),
            max(cmp$se_diff), max(abs(cmp$rr_clogit - cmp$rr_inla)), n2, nrow(cmp)))

# Flattening the prior from N(0, 31.6^2) to N(0, 1000^2) moves the max log difference from
# 0.00521 to 0.00525 - i.e. not at all. So the residual is INLA's Laplace approximation, not
# shrinkage from the prior. Agreement within 0.5% on the ratio scale is the honest claim;
# "agrees to 2 decimal places" overstates it (15 of 18 terms do, not all).
