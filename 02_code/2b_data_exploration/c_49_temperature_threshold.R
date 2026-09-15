############ A 4 degC THRESHOLD, FROM VIROLOGY RATHER THAN FROM WAIC ############
# HPAI survives for weeks in water and months in soil at or below about 4 degC, is preserved
# indefinitely when frozen, and is inactivated within days above roughly 20 degC. So if
# environmental persistence is the mechanism, the exposure that matters is not "temperature"
# but "how far below the survival threshold it was", and a threshold basis says that directly
# instead of asking a spline to discover it.
#
# THIS IS A DISCRIMINATING TEST, WHICH IS WHY IT IS WORTH RUNNING.
# c_46-c_48 found a WARM 90-day shock at lag 8-14 raising spillover where waterfowl had been
# abundant. Persistence virology predicts the opposite sign: COLD should help the virus survive.
# The two are different mechanisms and they make opposite predictions:
#
#   persistence    cold absolute conditions preserve virus in water and soil   -> cold harms
#   mobilisation   a late-winter thaw releases virus preserved in ice and snow,
#                  drives runoff, and coincides with waterfowl arrival         -> warm shock harms
#
# The shock result supports mobilisation. This tests persistence on its own terms: ABSOLUTE
# temperature, threshold at 4 degC, effect accruing as it gets colder. Both can be true - a cold
# winter banks the virus and a thaw releases it - but they are separate claims and only the
# absolute-scale threshold can speak to the first.
#
# Thresholds are fitted at 0, 4 and 10 degC. 10 is what the existing pipeline selected on WAIC
# for main effects (INLA_TEMP_THRESHOLD_C); 4 is the virological one. Reporting all three keeps
# it honest about how much the choice matters.

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
OTHER <- c("precipitation", "soil_moisture", "wind_speed", "runoff")

dat <- setDT(readRDS(paste0(objects_folder, "case_crossover_df_shock_post7.RDS")))
if (anyDuplicated(names(dat))) dat <- dat[, which(!duplicated(names(dat))), with = FALSE]
bl <- paste0("anseriformes_Lag", BIRDLAG)
dat[, bird_dep := rowMeans(.SD, na.rm = TRUE), .SDcols = bl]
dat[, birdhigh := as.integer(bird_dep > median(bird_dep, na.rm = TRUE)), by = stratum_id]

mask <- function(cb, keep) { m <- unclass(cb) * keep; attributes(m) <- attributes(cb)
                             dimnames(m) <- dimnames(cb); m }
lagw <- function(cb, lags, at, cen) {
  bs <- function(x, a) do.call(dlnm::onebasis, c(list(x = x), a))
  ev <- as.numeric(bs(at, attr(cb, "argvar")) - bs(cen, attr(cb, "argvar")))
  as.vector(t(outer(ev, colSums(bs(lags, attr(cb, "arglag")), na.rm = TRUE))))
}

run <- function(d, thr_C, panel) {
  # each lag column was z-scaled on its own mean/sd, so the Celsius cut has to be mapped per panel
  tz <- temp_threshold_z(d, thr_C)
  av <- list(fun = "thr", thr.value = tz, side = "l")   # effect accrues BELOW the threshold
  bases <- list(); fr <- data.frame(outbreak_binary = d$outbreak_binary)
  add <- function(nm, m) { bases[[nm]] <<- m; fr <<- cbind(fr, as.data.frame(unclass(m))) }

  cbT <- build_met_crossbasis(d, "temperature", MAXLAG, av, list(fun = "ns", df = 2))
  for (grp in c("high", "low")) {
    keep <- if (grp == "high") d$birdhigh else 1 - d$birdhigh
    m <- mask(cbT, keep); colnames(m) <- paste0("cb_temperature_", grp, ".", colnames(cbT))
    add(paste0("temperature_", grp), m)
  }
  for (v in OTHER) {
    cb <- build_met_crossbasis(d, v, MAXLAG, list(fun = "ns", df = 2), list(fun = "ns", df = 2))
    colnames(cb) <- paste0("cb_", v, ".", colnames(cb)); add(v, cb)
  }
  cbb <- build_met_crossbasis(d, "anseriformes", 28, list(fun = "lin"), list(fun = "ns", df = 2))
  colnames(cbb) <- paste0("cb_anseriformes.", colnames(cbb)); add("anseriformes", cbb)
  fr$id_stratum <- as.integer(as.factor(d$stratum_id))

  obj <- list(data = fr, bases = bases, cb_cols = lapply(bases, colnames), met_vars = names(bases),
              contrast_steps = list(), weekly = list(), weekly_vars = character(0),
              max_lag = MAXLAG, max_lag_by_var = list(), argvar = av, argvar_by_var = list(),
              arglag = list(fun = "ns", df = 2), arglag_by_var = list())
  fit <- suppressWarnings(fit_cc_inla_dlnm(obj, fixed_prec = PREC))

  # contrast: 5 degC COLDER than the threshold, against the threshold itself
  step_z <- 5 / sd(d$temperature, na.rm = TRUE)
  rbindlist(lapply(names(WINDOW), function(win) {
    lags <- WINDOW[[win]]
    est <- function(grp) {
      cb <- bases[[paste0("temperature_", grp)]]
      wt <- lagw(cb, lags, at = tz - step_z, cen = tz)
      cv <- dlnm_coef_vcov(fit, colnames(cb))
      lr <- sum(wt * cv$coef); se <- sqrt(as.numeric(t(wt) %*% cv$vcov %*% wt))
      c(lr = lr, se = se)
    }
    a <- est("high"); b <- est("low")
    cols <- c(obj$cb_cols[["temperature_high"]], obj$cb_cols[["temperature_low"]])
    cv <- dlnm_coef_vcov(fit, cols)
    w1 <- lagw(bases[["temperature_high"]], lags, at = tz - step_z, cen = tz)
    ct <- c(w1, -w1); lr <- sum(ct * cv$coef)
    se <- sqrt(as.numeric(t(ct) %*% cv$vcov %*% ct))
    f <- function(e) sprintf("%.2f (%.2f, %.2f)", exp(e[["lr"]]),
                             exp(e[["lr"]] - 1.96 * e[["se"]]), exp(e[["lr"]] + 1.96 * e[["se"]]))
    data.table(Panel = panel, Threshold = sprintf("%d degC", thr_C), Window = win,
               pct_below = sprintf("%.0f%%", 100 * mean(d$temperature < thr_C + 273.15, na.rm = TRUE)),
               WAIC = round(fit$waic$waic, 1),
               high = f(a), low = f(b),
               ratio = sprintf("%.2f (%.2f, %.2f)", exp(lr), exp(lr - 1.96 * se), exp(lr + 1.96 * se)),
               p_int = round(2 * pnorm(-abs(lr / se)), 4))
  }))
}

res <- rbindlist(lapply(c("Minnesota", "Pooled"), function(pn) {
  d <- if (pn == "Minnesota") dat[state == "Minnesota"] else copy(dat)
  rbindlist(lapply(c(0, 4, 10), function(th) {
    cat(sprintf("%-10s threshold %2d degC\n", pn, th)); flush.console()
    r <- try(run(d, th, pn), silent = TRUE)
    if (inherits(r, "try-error")) { cat("   failed:", attr(r, "condition")$message, "\n"); NULL } else r
  }))
}))
saveRDS(res, paste0(objects_folder, "si_temperature_threshold.RDS"))
cat("\n=== rate ratio per 5 degC BELOW the threshold, by waterfowl three weeks earlier ===\n")
cat("    (persistence predicts high > 1; the warm-shock result predicts the opposite)\n\n")
print(res[, .(Panel, Threshold, Window, pct_below, WAIC, high, low, ratio, p_int)], nrows = 60)
