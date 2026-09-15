############ SEASON INTERACTION TEST ############
# Is any exposure effect stronger in the spring wave than outside it?
#
# THE IDEA
# Season is the same for every day in a stratum, so it cannot be a plain covariate - the
# stratum intercept eats it. It can only act as a modifier. So we split each exposure into
# two columns:
#
#     x_spring  = x  in spring strata, 0 elsewhere
#     x_outside = x  outside spring,   0 in spring
#
# Put both in one model and each gets its own effect. The difference between them IS the
# interaction, and because both come from the same fit we can get a proper interval for it.
#
# WHY WE SPLIT THE CROSSBASIS AND NOT THE COLUMN
# Splitting the raw exposure column (x * spring) would be simpler, but it only works for a
# linear term: ns(x * 0) is not 0, so zeroing the input does not zero a spline basis. The
# primary model uses ns(df = 2) for meteorology, so we build the crossbasis first and mask
# THAT instead. Forcing everything linear to keep the code short changes the answer - it
# turned up 25 "significant" interactions against 7, almost all meteorological.

rm(list = ls())
project.folder = paste0(print(here::here()), '/')
source(paste0(project.folder, 'create_folder_structure.R'))
source(paste0(functions.folder, 'script_initiate.R'))
suppressMessages({library(INLA); library(Matrix)})
INLA::inla.setOption(num.threads = 1)
source(paste0(functions.folder, 'inla_dlnm_helpers.R'))

PREC   <- 4
LAGS   <- 0:28
WINDOW <- list(`0-7 days` = 0:7, `8-14 days` = 8:14,
               `15-21 days` = 15:21, `22-28 days` = 22:28)
# exposures come from INLA_PRIMARY_* in inla_dlnm_helpers.R - never redefined locally
VARS   <- INLA_PRIMARY_MET
LABEL <- INLA_PRIMARY_LAB

# ---- 1. panel, plus a spring flag that is constant within each stratum ----
dat <- setDT(readRDS(paste0(objects_folder, "case_crossover_df_shock_post7.RDS")))
if (anyDuplicated(names(dat))) dat <- dat[, which(!duplicated(names(dat))), with = FALSE]
dat[, date := as.Date(date)]
dat[, spring := as.integer(as.integer(strftime(date, "%m"))[outbreak_binary == 1][1] %in% 3:5),
    by = stratum_id]

# ---- 2. climatological anomaly columns (Minnesota only - that is all the 25-yr export has) ----
clim <- setDT(readRDS(paste0(objects_folder, "climate_anomaly_panel.RDS")))
clim <- clim[, c("zone_id", "week_start", paste0(VARS, "_z_25yr")), with = FALSE]
clim[, week_start := as.Date(week_start)]
setkey(clim, zone_id, week_start)
for (L in LAGS) {
  matched <- clim[.(dat$zone_id, dat$date - L), roll = TRUE, on = .(zone_id, week_start)]
  for (v in VARS) set(dat, j = sprintf("%s_clim_Lag%d", v, L),
                      value = as.numeric(scale(matched[[paste0(v, "_z_25yr")]])))
}

# ---- 3. one crossbasis, masked to a season ----
# Zero every row that belongs to the other season, keeping all the dlnm attributes so the
# contrast machinery downstream still knows what basis it is looking at.
mask_to_season <- function(cb, keep_row) {
  masked <- unclass(cb) * keep_row
  attributes(masked) <- attributes(cb)
  dimnames(masked)   <- dimnames(cb)
  masked
}

# The weights that turn crossbasis coefficients into "OR per +0.5 SD summed over these lags".
# It has to come from the basis the model was FITTED with - rebuilding a spline on each
# window's own range gives a different basis and a wrong (and window-invariant) answer.
window_weights <- function(obj, v, lags) {
  cb <- obj$bases[[v]]
  basis <- function(x, args) do.call(dlnm::onebasis, c(list(x = x), args))
  exposure_part <- as.numeric(basis(0.5, attr(cb, "argvar")) - basis(0, attr(cb, "argvar")))
  lag_part      <- colSums(basis(lags, attr(cb, "arglag")), na.rm = TRUE)
  as.vector(t(outer(exposure_part, lag_part)))
}

# ---- 4. fit one parameterisation on one panel, and read off the three quantities ----
run <- function(d, exposures, param, panel) {
  # exposure-response: ns(2) for meteorology, linear for waterfowl - same as the primary
  shape_of <- function(v) if (v == "anseriformes") list(fun = "lin") else list(fun = "ns", df = 2)
  lag_shape <- list(fun = "ns", df = 2)

  # build each variable's crossbasis once, then keep a spring copy and an outside copy
  bases <- list(); frame <- data.frame(outbreak_binary = d$outbreak_binary)
  for (v in exposures) {
    cb <- build_met_crossbasis(d, v, max(LAGS), shape_of(v), lag_shape)
    for (season in c("spring", "outside")) {
      keep <- if (season == "spring") d$spring else 1 - d$spring
      m <- mask_to_season(cb, keep)
      colnames(m) <- paste0("cb_", v, "_", season, ".", colnames(cb))
      bases[[paste0(v, "_", season)]] <- m
      frame <- cbind(frame, as.data.frame(unclass(m)))
    }
  }
  frame$id_stratum <- as.integer(as.factor(d$stratum_id))

  # the helper's fitting path expects this shape; we assemble it directly because
  # assemble_cc_dlnm has no notion of a masked basis
  obj <- list(data = frame, bases = bases,
              cb_cols = lapply(bases, colnames), met_vars = names(bases),
              contrast_steps = list(), weekly = list(), weekly_vars = character(0),
              max_lag = max(LAGS), max_lag_by_var = list(),
              argvar = list(fun = "ns", df = 2), argvar_by_var = list(),
              arglag = lag_shape, arglag_by_var = list())
  fit <- suppressWarnings(fit_cc_inla_dlnm(obj, fixed_prec = PREC))

  rbindlist(lapply(exposures, function(v) rbindlist(lapply(names(WINDOW), function(win) {
    lags <- WINDOW[[win]]
    a <- met_lag_effect(fit, obj, paste0(v, "_spring"),  lags)   # spring effect
    b <- met_lag_effect(fit, obj, paste0(v, "_outside"), lags)   # outside-spring effect

    # the interaction is the ratio of the two. on the log scale that is a difference, and
    # since the two terms are in one fit their covariance gives it a proper interval.
    cols <- c(obj$cb_cols[[paste0(v, "_spring")]], obj$cb_cols[[paste0(v, "_outside")]])
    cv   <- dlnm_coef_vcov(fit, cols)
    wt   <- window_weights(obj, paste0(v, "_spring"), lags)
    contrast <- c(wt, -wt)                    # +1 on the spring block, -1 on the other
    logratio <- sum(contrast * cv$coef)
    se       <- sqrt(as.numeric(t(contrast) %*% cv$vcov %*% contrast))

    data.table(parameterisation = param, panel = panel,
               exposure = unname(LABEL[sub("_shock$|_clim$", "", v)]), window = win,
               spring  = sprintf("%.2f (%.2f, %.2f)", a[["OR"]], a[["low"]], a[["high"]]),
               outside = sprintf("%.2f (%.2f, %.2f)", b[["OR"]], b[["low"]], b[["high"]]),
               ratio   = sprintf("%.2f (%.2f, %.2f)", exp(logratio),
                                 exp(logratio - 1.96 * se), exp(logratio + 1.96 * se)),
               p_interaction = round(2 * pnorm(-abs(logratio / se)), 3))
  }))))
}

# ---- 5. three exposure parameterisations ----
SETUPS <- list(
  list(param = "Absolute conditions",
       exposures = c(VARS, "anseriformes"),               panels = c("Minnesota", "Pooled")),
  list(param = "90-day rolling anomaly",
       exposures = c(paste0(VARS, "_shock"), "anseriformes"), panels = c("Minnesota", "Pooled")),
  list(param = "25-year climatological anomaly",
       exposures = c(paste0(VARS, "_clim"), "anseriformes"),  panels = "Minnesota"))

results <- rbindlist(lapply(SETUPS, function(s) {
  rbindlist(lapply(s$panels, function(p) {
    d <- copy(if (p == "Minnesota") dat[state == "Minnesota"] else dat)
    cat(sprintf("%-32s %-10s n = %d\n", s$param, p, sum(d$outbreak_binary)))
    run(d, s$exposures, s$param, p)
  }))
}))
saveRDS(results, paste0(objects_folder, "si_season_interaction.RDS"))
fwrite(results, file.path(tables_main_folder, "TableS7_season_interaction.csv"))

cat("\n=== interactions significant at 0.05 ===\n")
print(results[p_interaction < 0.05,
              .(parameterisation, panel, exposure, window, ratio, p_interaction)])
cat(sprintf("\n%d of %d tests significant (expect %.1f by chance)\n",
            sum(results$p_interaction < 0.05), nrow(results),
            0.05 * nrow(results)))
