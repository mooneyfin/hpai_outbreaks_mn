############ SPECIFICATION SENSITIVITY, ALL EXPOSURES AND ALL LAG WINDOWS ############
# The cached object stored one lag window per exposure, which forced the table to show
# meteorology at 0-7 d and waterfowl at 22-28 d - an arbitrary-looking choice, and one that hid
# whether a specification changes the lag profile rather than just the headline number.
#
# Refit keeping every exposure at every window, so the table can show the whole profile and
# nothing has to be picked. Two blocks: lag structure and prior. dWAIC belongs to the
# specification, so it is carried on the specification row rather than repeated per exposure.

rm(list = ls())
project.folder = paste0(print(here::here()), '/')
source(paste0(project.folder, 'create_folder_structure.R'))
source(paste0(functions.folder, 'script_initiate.R'))
suppressMessages({library(INLA); library(Matrix)})
INLA::inla.setOption(num.threads = 2)
source(paste0(functions.folder, 'inla_dlnm_helpers.R'))

NS2 <- list(fun = "ns", df = 2); LIN <- list(fun = "lin")
# exposures come from INLA_PRIMARY_* in inla_dlnm_helpers.R - never redefined locally
SHOCK <- INLA_PRIMARY_SHOCK
LEVEL <- INLA_PRIMARY_MET
MODELS <- list(`Absolute conditions` = LEVEL, `90-day shock` = SHOCK)
LAB <- INLA_PRIMARY_LAB
cell <- function(e) sprintf("%.2f (%.2f, %.2f)", e[["OR"]], e[["low"]], e[["high"]])

D <- setDT(readRDS(paste0(objects_folder, "case_crossover_df_shock_post7.RDS")))
if (anyDuplicated(names(D))) D <- D[, which(!duplicated(names(D))), with = FALSE]

# all three panels, not just pooled. The table used to report one panel without saying which,
# which left a reader unable to tell whether a sensitivity applied to Minnesota or everything.
PANELS <- c("Minnesota", "Other flyway states", "Pooled")
# DP is the panel currently being fitted. It has to live at top level: fit_spec's default
# argument resolves lexically from where the function was DEFINED, so a DP local to the loop
# below would be invisible to it.
DP <- NULL
# braces matter here: a top-level if/else split across lines parses as two statements and the
# else is orphaned
panel_of <- function(pn) {
  if (pn == "Minnesota") D[state == "Minnesota"]
  else if (pn == "Pooled") copy(D)
  else D[state != "Minnesota"]
}

fit_spec <- function(vars, max_lag = 28, lag_df = 2, prec = INLA_PREC_PRIMARY, D = DP) {
  mv <- c(vars, "anseriformes")
  av <- c(setNames(rep(list(NS2), length(vars)), vars), list(anseriformes = LIN))
  o <- assemble_cc_dlnm(D, met_vars = mv, weekly_vars = character(0), max_lag = max_lag,
                        argvar_by_var = av, arglag = list(fun = "ns", df = lag_df))
  f <- suppressWarnings(fit_cc_inla_dlnm(o, fixed_prec = prec))
  list(obj = o, fit = f, mv = mv, waic = f$waic$waic, max_lag = max_lag)
}
rows_of <- function(r, block, spec, mn, dwaic, diag = "") {
  rbindlist(lapply(r$mv, function(v) {
    lb <- unname(LAB[sub("_shock$", "", v)])
    vals <- vapply(names(INLA_LAG_WINDOWS), function(k) {
      w <- INLA_LAG_WINDOWS[[k]]
      if (w[2] > r$max_lag) return(NA_character_)      # window beyond this fit's lag range
      cell(met_lag_effect(r$fit, r$obj, v, seq(w[1], w[2])))
    }, character(1))
    c(list(Model = mn, Block = block, Specification = spec, dWAIC = dwaic,
           `Prior-implied OR range` = diag, Exposure = lb), as.list(vals))
  }))
}

LAGS <- list(`Lag 0-28 d, ns(df = 2)  [selected]` = list(m = 28, d = 2),
             `Lag 0-14 d, ns(df = 2)` = list(m = 14, d = 2),
             `Lag 0-21 d, ns(df = 2)` = list(m = 21, d = 2),
             `Lag 0-28 d, ns(df = 3)` = list(m = 28, d = 3),
             `Lag 0-28 d, ns(df = 4)` = list(m = 28, d = 4))
PRIORS <- c(0.001, 1, 4, 16, 64, 256)

set.seed(1)
out <- rbindlist(lapply(PANELS, function(pn) {
DP <<- panel_of(pn)
cat(sprintf("\n=== %s (%d cases) ===\n", pn, sum(DP$outbreak_binary)))
rbindlist(lapply(names(MODELS), function(mn) {
  vars <- MODELS[[mn]]
  ref <- fit_spec(vars)
  dw <- function(w) sprintf("%+.1f", w - ref$waic)

  lagrows <- rbindlist(lapply(names(LAGS), function(nm) {
    r <- fit_spec(vars, max_lag = LAGS[[nm]]$m, lag_df = LAGS[[nm]]$d)
    cat(sprintf("  %-20s lag  %-34s dWAIC %s\n", mn, nm, dw(r$waic)))
    rows_of(r, "Lag structure", nm, mn, dw(r$waic))
  }))

  # prior-implied OR: what the prior alone says about a +0.5 SD contrast, before any data
  cb <- ref$obj$bases[[vars[3]]]
  ob <- function(x, a) do.call(dlnm::onebasis, c(list(x = x), a))
  vv <- as.numeric(ob(0.5, attr(cb, "argvar")) - ob(0, attr(cb, "argvar")))
  cvec <- as.vector(t(outer(vv, colSums(ob(0:7, attr(cb, "arglag")), na.rm = TRUE))))

  priorrows <- rbindlist(lapply(PRIORS, function(pr) {
    r <- fit_spec(vars, prec = pr)
    OR <- exp(as.numeric(matrix(rnorm(20000 * length(cvec), 0, 1 / sqrt(pr)),
                                ncol = length(cvec)) %*% cvec))
    nm <- sprintf("N(0, %.3g^2)%s", 1 / sqrt(pr), if (pr == 4) "  [selected]" else "")
    cat(sprintf("  %-20s prior %-33s dWAIC %s\n", mn, nm, dw(r$waic)))
    rows_of(r, "Prior", nm, mn, dw(r$waic),
            sprintf("%.2f to %.2f", quantile(OR, .025), quantile(OR, .975)))
  }))
  cbind(Panel = pn, rbind(lagrows, priorrows))
}))
}))

saveRDS(out, paste0(objects_folder, "si_specification_full.RDS"))
cat("\nwrote si_specification_full.RDS |", nrow(out), "rows\n")
