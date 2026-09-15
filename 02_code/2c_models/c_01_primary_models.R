############ PRIMARY MODELS — SHOCK SPECIFICATION ############
# Final fits for the GeoHealth revision. Three panels (Minnesota, other flyway states,
# pooled) under the locked spec, plus the level-parameterisation companion that Figure 3a
# and SI S6 are built from. Everything downstream reads the cached fits from here.
#
# Spec and the reasoning behind each choice is in 07_drafts/REVISION_BUILD_PLAN.md.
# The one that needs restating here: the prior is N(0, 0.5^2), chosen on prior-predictive
# grounds (implied OR 0.66-1.51), NOT on WAIC. WAIC's own optimum sits at prec 64 and
# shrinks soil moisture to 0.91 and waterfowl to 1.79, which is shrinkage selection, not
# effect estimation.

rm(list = ls())
project.folder = paste0(print(here::here()), '/')
source(paste0(project.folder, 'create_folder_structure.R'))
source(paste0(functions.folder, 'script_initiate.R'))
suppressMessages({library(INLA); library(Matrix)})
INLA::inla.setOption(num.threads = 1)
source(paste0(functions.folder, 'inla_dlnm_helpers.R'))

PREC   <- INLA_PREC_PRIMARY                  # N(0, 0.5^2), defined once in inla_dlnm_helpers.R
NS2    <- list(fun = "ns", df = 2)
LIN    <- list(fun = "lin")
MAXLAG <- 28

# Runoff is IN, on all six exposures. It was out before because it is null as a main effect and
# costs WAIC (+7.6 pooled). It goes back in because the stated hypothesis names it and because a
# term that appears in the interaction model but not the primary is indefensible - the exposure
# set has to be the same object everywhere or the exhibits cannot be read against each other.
# The WAIC penalty is the price of that consistency, and it gets reported rather than avoided.
# exposures come from INLA_PRIMARY_* in inla_dlnm_helpers.R - never redefined locally
SHOCK_VARS   <- INLA_PRIMARY_SHOCK
LEVEL_VARS   <- INLA_PRIMARY_MET
BIRD         <- INLA_BIRD_VAR
SHOCK_LABELS <- INLA_PRIMARY_LAB
LEVEL_LABELS <- INLA_PRIMARY_LAB

# 1. Panels. The shock panel already carries the level columns, so both parameterisations
#    are fit on identical rows and their WAICs are directly comparable.
D <- setDT(readRDS(paste0(objects_folder, "case_crossover_df_shock_post7.RDS")))
if (anyDuplicated(names(D))) D <- D[, which(!duplicated(names(D))), with = FALSE]

PANELS <- list(
  Minnesota      = function(d) d[state == "Minnesota"],
  `Other flyway` = function(d) d[state != "Minnesota"],
  Pooled         = function(d) copy(d))

# 2. Fit one parameterisation on one panel.
fit_one <- function(dat, vars, labels) {
  av <- c(setNames(rep(list(NS2), length(vars)), vars), setNames(list(LIN), BIRD))
  o  <- assemble_cc_dlnm(dat, met_vars = c(vars, BIRD), weekly_vars = character(0),
                         max_lag = MAXLAG, argvar_by_var = av,
                         arglag = list(fun = "ns", df = 2))
  f  <- suppressWarnings(fit_cc_inla_dlnm(o, fixed_prec = PREC))
  list(obj = o, fit = f, labels = labels,
       cumulative  = lag_window_table(f, o, "cumulative",  labels = labels),
       lagspecific = lag_window_table(f, o, "lagspecific", labels = labels),
       prob        = prob_window_table(f, o, labels = labels),
       n = sum(dat$outbreak_binary), waic = f$waic$waic, pD = f$waic$p.eff)
}

fits <- list()
for (pn in names(PANELS)) {
  dp <- PANELS[[pn]](D)
  for (kind in c("shock", "level")) {
    vars <- if (kind == "shock") SHOCK_VARS else LEVEL_VARS
    labs <- if (kind == "shock") SHOCK_LABELS else LEVEL_LABELS
    cat(sprintf("fitting %-13s | %-5s | n = %d ... ", pn, kind, sum(dp$outbreak_binary)))
    r <- fit_one(dp, vars, labs)
    fits[[paste(pn, kind, sep = "|")]] <- r
    cat(sprintf("WAIC %.1f (pD %.1f)\n", r$waic, r$pD))
  }
}

saveRDS(fits, paste0(objects_folder, "primary_shock_level_fits.RDS"))

# 3. Print the headline tables so the log is self-documenting.
for (nm in names(fits)) {
  r <- fits[[nm]]
  cat(sprintf("\n================ %s | n = %d | WAIC %.1f ================\n", nm, r$n, r$waic))
  cat("-- cumulative --\n");   print(r$cumulative)
  cat("-- lag-specific --\n"); print(r$lagspecific)
}
