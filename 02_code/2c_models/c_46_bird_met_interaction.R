############ DOES METEOROLOGY ONLY MATTER WHERE WATERFOWL WERE? ############
# The hypothesised pathway: waterfowl deposit virus around three weeks out, then conditions in
# the fortnight before detection decide whether it survives long enough to get into a barn.
#
# THIS IS AN INTERACTION, NOT MEDIATION. A mediator has to be CAUSED BY the exposure, and weather
# is not caused by bird abundance. The real mediator here is environmental virus load, which is
# unmeasured, so no mediation estimator - g-computation included - can identify a natural indirect
# effect. What IS testable is whether the meteorological effect differs by how many birds were
# around three weeks earlier.
#
# WHY IT MATTERS THAT WE HAVE NOT TESTED IT. The primary model enters birds and meteorology as
# additive terms on the log scale, i.e. independent multiplicative effects. If weather only acts
# where virus was deposited, the primary is averaging the weather effect over cell-days where it
# CANNOT act and cell-days where it can, which dilutes it toward the null. That is a live
# mechanism for a false negative and it is the last one standing.
#
# THE SPLIT. Anseriformes averaged over lag 21-28, dichotomised at the median WITHIN each stratum.
# Within-stratum because a stratum is one cell in one calendar month, so an above-median day is
# "more birds about three weeks before this day than on other days of that month here" - which is
# the contrast the case-crossover identifies anyway. A panel-wide median would mostly encode
# season, and we already tested season (Table S6).

project.folder = paste0(print(here::here()), '/')
src <- readLines(paste0(project.folder, 'create_folder_structure.R'))
eval(parse(text = paste(src[!grepl("^\\s*rm\\(list", src)], collapse = "\n")))
source(paste0(functions.folder, 'script_initiate.R'))
suppressMessages({library(INLA); library(Matrix)})
INLA::inla.setOption(num.threads = 4)
source(paste0(functions.folder, 'inla_dlnm_helpers.R'))

PREC   <- 4
MAXLAG <- 14                      # the hypothesised weather window
WINDOW <- list(`0-7 days` = 0:7, `8-14 days` = 8:14)
BIRDLAG <- 21:28                  # the hypothesised deposition window
VARS   <- c("temperature", "precipitation", "soil_moisture", "wind_speed")
LAB    <- setNames(c("Temperature", "Precipitation", "Soil moisture", "Wind speed"), VARS)

dat <- setDT(readRDS(paste0(objects_folder, "case_crossover_df_shock_post7.RDS")))
if (anyDuplicated(names(dat))) dat <- dat[, which(!duplicated(names(dat))), with = FALSE]
dat[, date := as.Date(date)]

# birds three weeks out, then split within stratum
bl <- paste0("anseriformes_Lag", BIRDLAG)
dat[, bird_dep := rowMeans(.SD, na.rm = TRUE), .SDcols = bl]
dat[, birdhigh := as.integer(bird_dep > median(bird_dep, na.rm = TRUE)), by = stratum_id]
cat(sprintf("rows %d | high-bird rows %d (%.0f%%) | cases in high-bird rows %d of %d\n",
            nrow(dat), sum(dat$birdhigh), 100 * mean(dat$birdhigh),
            sum(dat$outbreak_binary[dat$birdhigh == 1]), sum(dat$outbreak_binary)))
# a stratum with no contrast contributes nothing to the interaction
dat[, nb := uniqueN(birdhigh), by = stratum_id]
cat(sprintf("strata with both high and low rows: %d of %d\n",
            uniqueN(dat[nb == 2]$stratum_id), uniqueN(dat$stratum_id)))

# zero the rows belonging to the other group, keeping dlnm's attributes intact
mask <- function(cb, keep) { m <- unclass(cb) * keep; attributes(m) <- attributes(cb)
                             dimnames(m) <- dimnames(cb); m }
window_weights <- function(obj, v, lags) {
  cb <- obj$bases[[v]]
  bs <- function(x, a) do.call(dlnm::onebasis, c(list(x = x), a))
  ev <- as.numeric(bs(0.5, attr(cb, "argvar")) - bs(0, attr(cb, "argvar")))
  as.vector(t(outer(ev, colSums(bs(lags, attr(cb, "arglag")), na.rm = TRUE))))
}

run <- function(d, mvars, panel, kind) {
  bases <- list(); frame <- data.frame(outbreak_binary = d$outbreak_binary)
  for (v in mvars) {
    cb <- build_met_crossbasis(d, v, MAXLAG, list(fun = "ns", df = 2), list(fun = "ns", df = 2))
    for (grp in c("high", "low")) {
      keep <- if (grp == "high") d$birdhigh else 1 - d$birdhigh
      m <- mask(cb, keep); colnames(m) <- paste0("cb_", v, "_", grp, ".", colnames(cb))
      bases[[paste0(v, "_", grp)]] <- m
      frame <- cbind(frame, as.data.frame(unclass(m)))
    }
  }
  # waterfowl stays in as its own term so the interaction is not just picking up the bird effect
  cbb <- build_met_crossbasis(d, "anseriformes", 28, list(fun = "lin"), list(fun = "ns", df = 2))
  colnames(cbb) <- paste0("cb_anseriformes.", colnames(cbb))
  bases[["anseriformes"]] <- cbb; frame <- cbind(frame, as.data.frame(unclass(cbb)))
  frame$id_stratum <- as.integer(as.factor(d$stratum_id))

  obj <- list(data = frame, bases = bases, cb_cols = lapply(bases, colnames),
              met_vars = names(bases), contrast_steps = list(), weekly = list(),
              weekly_vars = character(0), max_lag = MAXLAG, max_lag_by_var = list(),
              argvar = list(fun = "ns", df = 2), argvar_by_var = list(),
              arglag = list(fun = "ns", df = 2), arglag_by_var = list())
  fit <- suppressWarnings(fit_cc_inla_dlnm(obj, fixed_prec = PREC))

  rbindlist(lapply(mvars, function(v) rbindlist(lapply(names(WINDOW), function(win) {
    lags <- WINDOW[[win]]
    a <- met_lag_effect(fit, obj, paste0(v, "_high"), lags)
    b <- met_lag_effect(fit, obj, paste0(v, "_low"),  lags)
    cols <- c(obj$cb_cols[[paste0(v, "_high")]], obj$cb_cols[[paste0(v, "_low")]])
    cv <- dlnm_coef_vcov(fit, cols)
    wt <- window_weights(obj, paste0(v, "_high"), lags)
    ct <- c(wt, -wt)                       # high minus low, on the log scale
    lr <- sum(ct * cv$coef)
    se <- sqrt(as.numeric(t(ct) %*% cv$vcov %*% ct))
    data.table(Panel = panel, Exposure = unname(LAB[sub("_shock$", "", v)]),
               Window = win, Kind = kind, n = sum(d$outbreak_binary),
               high = sprintf("%.2f (%.2f, %.2f)", a[["OR"]], a[["low"]], a[["high"]]),
               low  = sprintf("%.2f (%.2f, %.2f)", b[["OR"]], b[["low"]], b[["high"]]),
               ratio = sprintf("%.2f (%.2f, %.2f)", exp(lr), exp(lr - 1.96 * se), exp(lr + 1.96 * se)),
               p_int = round(2 * pnorm(-abs(lr / se)), 3))
  }))))
}

SETUPS <- list(list(kind = "Absolute conditions", vars = VARS),
               list(kind = "90-day shock",        vars = paste0(VARS, "_shock")))
res <- rbindlist(lapply(SETUPS, function(s) rbindlist(lapply(c("Minnesota", "Pooled"), function(pn) {
  d <- if (pn == "Minnesota") dat[state == "Minnesota"] else copy(dat)
  cat(sprintf("%-22s %-10s n = %d\n", s$kind, pn, sum(d$outbreak_binary))); flush.console()
  run(d, s$vars, pn, s$kind)
}))))
saveRDS(res, paste0(objects_folder, "si_bird_met_interaction.RDS"))

cat("\n=== meteorological effect, high vs low waterfowl three weeks earlier ===\n")
print(res[, .(Kind, Panel, Exposure, Window, high, low, ratio, p_int)], nrows = 100)
cat(sprintf("\ninteractions significant at 0.05: %d of %d (expect %.1f by chance)\n",
            sum(res$p_int < 0.05), nrow(res), 0.05 * nrow(res)))
if (any(res$p_int < 0.05)) print(res[p_int < 0.05, .(Kind, Panel, Exposure, Window, ratio, p_int)])
