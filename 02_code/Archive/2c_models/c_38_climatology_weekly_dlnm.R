############ CLIMATOLOGICAL ANOMALY AS A WEEKLY DLNM, WITH THE BASIS RESELECTED ############
# The 25-year climatology is weekly - one mean and sd per cell per week_idx - so every day
# inside a week is compared against the SAME normal. A 29-knot daily lag basis is therefore
# estimating structure finer than the exposure has, and spends 10 parameters per exposure doing
# it on 175 cases. Binning the lag axis to weeks matches the exposure's own resolution and
# drops that to 4 or fewer.
#
# The crossbasis cannot be carried over from the daily model unchanged, so it is reselected here
# by WAIC over exposure-response shape, lag shape and window length.

rm(list = ls())
project.folder = paste0(print(here::here()), '/')
source(paste0(project.folder, 'create_folder_structure.R'))
source(paste0(functions.folder, 'script_initiate.R'))
suppressMessages({library(INLA); library(Matrix)})
INLA::inla.setOption(num.threads = 2)
source(paste0(functions.folder, 'inla_dlnm_helpers.R'))

PREC <- 4
VARS <- c("temperature", "precipitation", "soil_moisture", "wind_speed")
CLIM <- paste0(VARS, "_clim")
LAB  <- c(setNames(c("Temperature", "Precipitation", "Soil moisture", "Wind speed"), CLIM),
          anseriformes = "Anseriformes")
cell <- function(e) sprintf("%.2f (%.2f, %.2f)", e[["OR"]], e[["low"]], e[["high"]])

D <- setDT(readRDS(paste0(objects_folder, "case_crossover_df_climanom_post7.RDS")))
if (anyDuplicated(names(D))) D <- D[, which(!duplicated(names(D))), with = FALSE]
have <- CLIM[vapply(CLIM, function(v) paste0(v, "_Lag0") %in% names(D), TRUE)]
cat("exposures available:", paste(have, collapse = ", "), "\n")

# weekly means from the daily lag columns, z-scored the same way the daily ones were
week_mat <- function(dat, v, nweek) {
  sapply(seq_len(nweek) - 1L, function(w) {
    cols <- paste0(v, "_Lag", (w * 7):(w * 7 + 6))
    cols <- intersect(cols, names(dat))
    as.numeric(scale(rowMeans(dat[, ..cols], na.rm = TRUE)))
  })
}
fit_weekly <- function(dat, vars, nweek, av, al) {
  bases <- list(); dd <- data.frame(outbreak_binary = dat$outbreak_binary)
  for (v in c(vars, "anseriformes")) {
    a <- if (v == "anseriformes") list(fun = "lin") else av
    cb <- dlnm::crossbasis(week_mat(dat, v, nweek), lag = c(0, nweek - 1),
                           argvar = a, arglag = al)
    colnames(cb) <- paste0("cb_", v, ".", colnames(cb))
    bases[[v]] <- cb; dd <- cbind(dd, as.data.frame(unclass(cb)))
  }
  dd$id_stratum <- as.integer(as.factor(dat$stratum_id))
  obj <- list(data = dd, bases = bases, cb_cols = lapply(bases, colnames),
              met_vars = names(bases), contrast_steps = list(), weekly = list(),
              weekly_vars = character(0), max_lag = nweek - 1, max_lag_by_var = list(),
              argvar = av, argvar_by_var = list(), arglag = al, arglag_by_var = list())
  list(obj = obj, fit = suppressWarnings(fit_cc_inla_dlnm(obj, fixed_prec = PREC)))
}

# ---- basis selection by WAIC, pooled panel ----
GRID <- CJ(nweek = c(4L, 8L),
           var_shape = c("lin", "ns2"),
           lag_shape = c("lin", "ns2"), sorted = FALSE)
shape <- function(s) if (s == "lin") list(fun = "lin") else list(fun = "ns", df = 2)

sel <- rbindlist(lapply(seq_len(nrow(GRID)), function(i) {
  g <- GRID[i]
  # ns(2) on a 4-week lag axis has only 4 distinct positions; skip combinations that can't fit
  r <- try(fit_weekly(D, have, g$nweek, shape(g$var_shape), shape(g$lag_shape)), silent = TRUE)
  if (inherits(r, "try-error")) return(NULL)
  cat(sprintf("  lag 0-%d wk | argvar %-3s | arglag %-3s | WAIC %.1f\n",
              g$nweek - 1, g$var_shape, g$lag_shape, r$fit$waic$waic))
  data.table(g, WAIC = r$fit$waic$waic)
}))
setorder(sel, WAIC)
best <- sel[1]
cat(sprintf("\nselected: lag 0-%d weeks, argvar %s, arglag %s (WAIC %.1f)\n",
            best$nweek - 1, best$var_shape, best$lag_shape, best$WAIC))

# ---- fit the selected spec on all three panels ----
PANELS <- c("Minnesota", "Other flyway states", "Pooled")
pick <- function(pn) {
  if (pn == "Minnesota") D[state == "Minnesota"]
  else if (pn == "Pooled") copy(D) else D[state != "Minnesota"]
}
out <- rbindlist(lapply(PANELS, function(pn) {
  d <- pick(pn)
  r <- fit_weekly(d, have, best$nweek, shape(best$var_shape), shape(best$lag_shape))
  rbindlist(lapply(names(r$obj$bases), function(v)
    rbindlist(lapply(seq_len(best$nweek) - 1L, function(w) data.table(
      Panel = pn, n = sum(d$outbreak_binary), Exposure = unname(LAB[v]),
      `Lag week` = paste0("Week ", w),
      est = cell(met_lag_effect(r$fit, r$obj, v, w)))))))
}))

saveRDS(list(selection = sel, estimates = out),
        paste0(objects_folder, "si_climatology_weekly.RDS"))
cat("\n=== selected weekly climatological DLNM ===\n")
print(dcast(out, Exposure + Panel ~ `Lag week`, value.var = "est"))
