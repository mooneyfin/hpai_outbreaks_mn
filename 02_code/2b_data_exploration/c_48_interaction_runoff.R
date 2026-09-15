############ THE HYPOTHESISED TRIPLE, WITH RUNOFF PUT BACK ############
# c_46/c_47 established that the meteorological effect at lag 8-14 days appears only where
# waterfowl were abundant three weeks earlier, and that it survives an unsmoothed lag basis and
# dropping the dichotomisation entirely. Those models did not contain runoff - it was cut from the
# primary because it was null as a MAIN effect.
#
# But "null as a main effect" is exactly what temperature looked like too, and for the same reason:
# averaged over cell-days where there was no virus to keep alive. Runoff belongs in the stated
# mechanism (wind, temperature, runoff), so it gets tested the same way here.
#
# Two exposure sets: the hypothesised triple on its own, and the full set with runoff added, so it
# is clear whether any runoff signal needs the other water terms out of the way.
# Only the two specifications that c_47 showed were the conservative ones: an unsmoothed strata
# lag basis, and the continuous modifier with no dichotomisation.

project.folder = paste0(print(here::here()), '/')
src <- readLines(paste0(project.folder, 'create_folder_structure.R'))
eval(parse(text = paste(src[!grepl("^\\s*rm\\(list", src)], collapse = "\n")))
source(paste0(functions.folder, 'script_initiate.R'))
suppressMessages({library(INLA); library(Matrix)})
INLA::inla.setOption(num.threads = 4)
source(paste0(functions.folder, 'inla_dlnm_helpers.R'))

PREC <- 4; MAXLAG <- 14
WINDOW <- list(`0-7 days` = 0:7, `8-14 days` = 8:14)
BIRDLAG <- 21:28
TRIPLE <- c("temperature", "runoff", "wind_speed")
FULL   <- c("temperature", "precipitation", "soil_moisture", "wind_speed", "runoff")
VARS   <- FULL
LAB  <- setNames(c("Temperature", "Precipitation", "Soil moisture", "Wind speed", "Runoff"), FULL)
EXPSETS <- list(`Triple (temp, runoff, wind)` = TRIPLE, `Full set + runoff` = FULL)

dat <- setDT(readRDS(paste0(objects_folder, "case_crossover_df_shock_post7.RDS")))
if (anyDuplicated(names(dat))) dat <- dat[, which(!duplicated(names(dat))), with = FALSE]
bl <- paste0("anseriformes_Lag", BIRDLAG)
dat[, bird_dep := rowMeans(.SD, na.rm = TRUE), .SDcols = bl]
dat[, bird_c := bird_dep - mean(bird_dep, na.rm = TRUE), by = stratum_id]   # within-stratum
dat[, bird_z := bird_c / sd(bird_c, na.rm = TRUE)]
dat[, birdhigh := as.integer(bird_dep > median(bird_dep, na.rm = TRUE)), by = stratum_id]
# rank-based tertiles, not cut() on quantiles: small strata have tied bird values, which makes
# the quantile breaks non-unique and cut() throws
dat[, brank := data.table::frank(bird_dep, ties.method = "average") / .N, by = stratum_id]
dat[, btert := fifelse(brank <= 1/3, 1L, fifelse(brank > 2/3, 3L, 2L))]

mask <- function(cb, keep) { m <- unclass(cb) * keep; attributes(m) <- attributes(cb)
                             dimnames(m) <- dimnames(cb); m }
scale_cb <- function(cb, z) { m <- unclass(cb) * z; attributes(m) <- attributes(cb)
                              dimnames(m) <- dimnames(cb); m }
lagw <- function(cb, lags) {
  bs <- function(x, a) do.call(dlnm::onebasis, c(list(x = x), a))
  ev <- as.numeric(bs(0.5, attr(cb, "argvar")) - bs(0, attr(cb, "argvar")))
  as.vector(t(outer(ev, colSums(bs(lags, attr(cb, "arglag")), na.rm = TRUE))))
}
ARGLAG <- list(`strata at 8` = list(fun = "strata", breaks = 8))

# assemble whichever design is asked for and read off the contrast
run <- function(d, mvars, arglag, mode, tag, panel, kind) {
  bases <- list(); fr <- data.frame(outbreak_binary = d$outbreak_binary)
  add <- function(nm, m) { bases[[nm]] <<- m; fr <<- cbind(fr, as.data.frame(unclass(m))) }
  for (v in mvars) {
    cb <- build_met_crossbasis(d, v, MAXLAG, list(fun = "ns", df = 2), arglag)
    if (mode == "continuous") {
      m1 <- cb; colnames(m1) <- paste0("cb_", v, "_main.", colnames(cb)); add(paste0(v, "_main"), m1)
      m2 <- scale_cb(cb, d$bird_z); colnames(m2) <- paste0("cb_", v, "_int.", colnames(cb))
      add(paste0(v, "_int"), m2)
    } else {
      for (grp in c("high", "low")) {
        keep <- if (grp == "high") d$grp_hi else 1 - d$grp_hi
        m <- mask(cb, keep); colnames(m) <- paste0("cb_", v, "_", grp, ".", colnames(cb))
        add(paste0(v, "_", grp), m)
      }
    }
  }
  cbb <- build_met_crossbasis(d, "anseriformes", 28, list(fun = "lin"), list(fun = "ns", df = 2))
  colnames(cbb) <- paste0("cb_anseriformes.", colnames(cbb)); add("anseriformes", cbb)
  if (mode == "continuous") fr$bird_z <- d$bird_z
  fr$id_stratum <- as.integer(as.factor(d$stratum_id))
  obj <- list(data = fr, bases = bases, cb_cols = lapply(bases, colnames), met_vars = names(bases),
              contrast_steps = list(), weekly = list(), weekly_vars = character(0),
              max_lag = MAXLAG, max_lag_by_var = list(), argvar = list(fun = "ns", df = 2),
              argvar_by_var = list(), arglag = arglag, arglag_by_var = list())
  fit <- suppressWarnings(fit_cc_inla_dlnm(obj, fixed_prec = PREC))

  rbindlist(lapply(mvars, function(v) rbindlist(lapply(names(WINDOW), function(win) {
    lags <- WINDOW[[win]]
    if (mode == "continuous") {
      cb <- bases[[paste0(v, "_int")]]; wt <- lagw(cb, lags)
      cv <- dlnm_coef_vcov(fit, colnames(cb))
      # ratio of the effect at +1 SD birds to the effect at -1 SD = exp(2 * w'b_int)
      lr <- 2 * sum(wt * cv$coef)
      se <- 2 * sqrt(as.numeric(t(wt) %*% cv$vcov %*% wt))
      hi <- lo <- ""
    } else {
      a <- met_lag_effect(fit, obj, paste0(v, "_high"), lags)
      b <- met_lag_effect(fit, obj, paste0(v, "_low"),  lags)
      hi <- sprintf("%.2f (%.2f, %.2f)", a[["OR"]], a[["low"]], a[["high"]])
      lo <- sprintf("%.2f (%.2f, %.2f)", b[["OR"]], b[["low"]], b[["high"]])
      cols <- c(obj$cb_cols[[paste0(v, "_high")]], obj$cb_cols[[paste0(v, "_low")]])
      cv <- dlnm_coef_vcov(fit, cols); w1 <- lagw(bases[[paste0(v, "_high")]], lags)
      ct <- c(w1, -w1); lr <- sum(ct * cv$coef)
      se <- sqrt(as.numeric(t(ct) %*% cv$vcov %*% ct))
    }
    data.table(Check = tag, Panel = panel, Kind = kind,
               Exposure = unname(LAB[sub("_shock$", "", v)]), Window = win,
               high = hi, low = lo,
               ratio = sprintf("%.2f (%.2f, %.2f)", exp(lr), exp(lr - 1.96 * se), exp(lr + 1.96 * se)),
               p_int = round(2 * pnorm(-abs(lr / se)), 4))
  }))))
}

grid <- CJ(set = names(EXPSETS), kind = c("Absolute conditions", "90-day shock"),
           panel = c("Minnesota", "Pooled"), sorted = FALSE)
res <- rbindlist(apply(grid, 1, function(g) {
  setnm <- g[["set"]]; kind <- g[["kind"]]; panel <- g[["panel"]]
  base <- EXPSETS[[setnm]]
  mvars <- if (kind == "90-day shock") paste0(base, "_shock") else base
  d0 <- if (panel == "Minnesota") dat[state == "Minnesota"] else copy(dat)
  out <- list()
  for (nm in c("Median split", "Continuous (no split)")) {
    d <- copy(d0)[, grp_hi := birdhigh]
    mode <- if (nm == "Continuous (no split)") "continuous" else "split"
    cat(sprintf("%-28s %-20s %-10s %-22s\n", setnm, kind, panel, nm)); flush.console()
    r <- try(run(d, mvars, ARGLAG[["strata at 8"]], mode, nm, panel, kind), silent = TRUE)
    if (inherits(r, "try-error")) { cat("   failed\n"); next }
    out[[nm]] <- cbind(Set = setnm, r)
  }
  rbindlist(out)
}))
saveRDS(res, paste0(objects_folder, "si_interaction_runoff.RDS"))

cat("\n================ lag 8-14 days, all exposures ================\n")
print(res[Window == "8-14 days", .(Set, Panel, Kind, Check, Exposure, high, low, ratio, p_int)],
      nrows = 200)
cat("\n================ runoff only ================\n")
print(res[Exposure == "Runoff", .(Set, Panel, Kind, Check, Window, ratio, p_int)], nrows = 100)
