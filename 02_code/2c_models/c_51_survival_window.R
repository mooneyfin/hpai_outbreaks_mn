############ WATERFOWL EFFECT MODIFIED BY TEMPERATURE IN THE SURVIVAL WINDOW ############
# Kurmi et al. (2013, Indian J Virol 24:272-77) measured H5N1 survival in poultry faeces:
# 8 weeks at 4 degC, 5 days at 24 degC, 24 h at 37 degC, 18 h at 42 degC. Taking logs, that is
# about -0.12 log-hours of survival per degC and it is near-constant from 4 to 37 degC.
#
# TWO THINGS FOLLOW, AND BOTH CHANGE THE MODEL.
#
# 1. NO THRESHOLD. 4 degC is not a special point, it is where they happened to measure. Survival
#    is log-linear in temperature across the range, so a LINEAR temperature term is the
#    mechanistically correct shape - simpler than the spline in c_46-c_48 and the threshold in
#    c_49, and much easier to defend than either.
#
# 2. WE HAD THE CONTRAST BACKWARDS. Their mechanism is: birds deposit virus, then the temperature
#    BETWEEN deposition and detection decides whether it is still viable. So the exposure is
#    waterfowl at lag 21-28 and the modifier is temperature over lag 0-21 - the survival window -
#    not meteorology at lag 8-14 modified by birds. An interaction is algebraically symmetric, so
#    it is the same product term, but the contrast reported and the window used were both wrong
#    for the biology.
#
# PREDICTION, STATED BEFORE FITTING: the waterfowl effect is LARGER when the intervening three
# weeks were cold. At 4 degC virus deposited three weeks ago is still viable; at 24 degC it died
# two weeks ago.
#
# This does not contradict the warm-shock result in c_50. A cold spell preserves, a thaw
# mobilises; they are different halves of the same story and they live in different windows.
#
# Their other finding - dry versus wet faeces barely differed at any temperature - says
# persistence is a temperature story, not a moisture one. So the protective precipitation effect
# is probably washout or dilution and should not be written up as persistence.

project.folder = paste0(print(here::here()), '/')
src <- readLines(paste0(project.folder, 'create_folder_structure.R'))
eval(parse(text = paste(src[!grepl("^\\s*rm\\(list", src)], collapse = "\n")))
source(paste0(functions.folder, 'script_initiate.R'))
suppressMessages({library(INLA); library(Matrix)})
INLA::inla.setOption(num.threads = 4)
source(paste0(functions.folder, 'inla_dlnm_helpers.R'))

PREC     <- INLA_PREC_PRIMARY
SURV_WIN <- 0:21          # deposition-to-detection window
BIRDWIN  <- 21:28         # deposition window
NSAMP    <- 4000
MET      <- c("precipitation", "soil_moisture", "wind_speed", "runoff")

cc <- setDT(readRDS(paste0(objects_folder, "case_crossover_df_shock_post7.RDS")))
if (anyDuplicated(names(cc))) cc <- cc[, which(!duplicated(names(cc))), with = FALSE]
cc[, date := as.Date(date)]

# The lag columns were each z-scaled on their own mean/sd, so they cannot give a temperature in
# degrees. Rebuild the survival-window mean from the raw daily panel instead.
daily <- setDT(readRDS(paste0(objects_folder, multistate_ts_daily_rds)))
daily[, date := as.Date(date)]
setorder(daily, zone_id, date)
daily[, tC := temperature - 273.15]
daily[, t_surv := data.table::frollmean(tC, length(SURV_WIN), align = "right"), by = zone_id]
cc <- merge(cc, daily[, .(zone_id, date, t_surv)], by = c("zone_id", "date"), all.x = TRUE)
cc <- cc[!is.na(t_surv)]
cat(sprintf("panel: %d cases | survival-window temp %.1f to %.1f degC (median %.1f)\n",
            sum(cc$outbreak_binary), min(cc$t_surv), max(cc$t_surv), median(cc$t_surv)))

# two versions of the modifier: absolute (what the virology is about) and within-stratum centred
cc[, t_abs := t_surv]
cc[, t_wsc := t_surv - mean(t_surv), by = stratum_id]
QS <- quantile(cc$t_abs, c(0.1, 0.9))
cat(sprintf("modifier contrast: %.1f degC (p10) vs %.1f degC (p90)\n", QS[1], QS[2]))

scale_cb <- function(cb, z) { m <- unclass(cb) * z; attributes(m) <- attributes(cb)
                              dimnames(m) <- dimnames(cb); m }
lagw <- function(cb, lags) {
  bs <- function(x, a) do.call(dlnm::onebasis, c(list(x = x), a))
  ev <- as.numeric(bs(0.5, attr(cb, "argvar")) - bs(0, attr(cb, "argvar")))
  as.vector(t(outer(ev, colSums(bs(lags, attr(cb, "arglag")), na.rm = TRUE))))
}

run <- function(d, panel, modkind, kind) {
  z <- if (modkind == "Absolute") d$t_abs else d$t_wsc
  sdz <- sd(z); zc <- (z - mean(z)) / sdz            # per +1 SD of the modifier
  mvars <- if (kind == "90-day shock") paste0(MET, "_shock") else MET
  bases <- list(); fr <- data.frame(outbreak_binary = d$outbreak_binary)
  add <- function(nm, m) { bases[[nm]] <<- m; fr <<- cbind(fr, as.data.frame(unclass(m))) }

  # waterfowl, and waterfowl x temperature. linear on the exposure scale, as in the primary.
  cbb <- build_met_crossbasis(d, "anseriformes", 28, list(fun = "lin"), list(fun = "ns", df = 2))
  m1 <- cbb; colnames(m1) <- paste0("cb_bird_main.", colnames(cbb)); add("bird_main", m1)
  m2 <- scale_cb(cbb, zc); colnames(m2) <- paste0("cb_bird_int.", colnames(cbb)); add("bird_int", m2)
  # temperature itself stays in the model; linear, per Kurmi's log-linear survival law
  cbT <- build_met_crossbasis(d, if (kind == "90-day shock") "temperature_shock" else "temperature",
                              21, list(fun = "lin"), list(fun = "ns", df = 2))
  colnames(cbT) <- paste0("cb_temp.", colnames(cbT)); add("temp", cbT)
  for (v in mvars) {
    cb <- build_met_crossbasis(d, v, 14, list(fun = "ns", df = 2), list(fun = "ns", df = 2))
    colnames(cb) <- paste0("cb_", v, ".", colnames(cb)); add(v, cb)
  }
  fr$tmod <- zc
  fr$id_stratum <- as.integer(as.factor(d$stratum_id))
  obj <- list(data = fr, bases = bases, cb_cols = lapply(bases, colnames), met_vars = names(bases),
              contrast_steps = list(), weekly = list(), weekly_vars = character(0),
              max_lag = 28, max_lag_by_var = list(), argvar = list(fun = "lin"),
              argvar_by_var = list(), arglag = list(fun = "ns", df = 2), arglag_by_var = list())
  fit <- suppressWarnings(fit_cc_inla_dlnm(obj, fixed_prec = PREC))

  samp <- INLA::inla.posterior.sample(NSAMP, fit, selection = list(), verbose = FALSE)
  nm_all <- sub(":[0-9]+$", "", rownames(samp[[1]]$latent))
  draws <- function(cols) t(vapply(samp, function(s) s$latent[match(cols, nm_all), 1],
                                   numeric(length(cols))))
  summ <- function(x) { q <- quantile(x, c(.5, .025, .975))
    list(est = sprintf("%.2f (%.2f, %.2f)", exp(q[1]), exp(q[2]), exp(q[3])),
         p = round(max(mean(x > 0), mean(x < 0)), 3)) }

  wt <- lagw(bases[["bird_main"]], BIRDWIN)
  dm <- draws(colnames(bases[["bird_main"]])) %*% wt
  di <- draws(colnames(bases[["bird_int"]]))  %*% wt
  # evaluate at the 10th and 90th percentile of the modifier, in SD units
  z10 <- (QS[1] - mean(z)) / sdz; z90 <- (QS[2] - mean(z)) / sdz
  cold <- summ(dm + di * z10); warm <- summ(dm + di * z90)
  rat  <- summ(di * (z10 - z90))                    # cold vs warm
  data.table(Panel = panel, Modifier = modkind, Kind = kind, n = sum(d$outbreak_binary),
             `Bird effect, cold` = cold$est, `P(dir)` = cold$p,
             `Bird effect, warm` = warm$est, `P(dir) ` = warm$p,
             `Cold vs warm` = rat$est, `P(dir)  ` = rat$p)
}

grid <- CJ(kind = c("Absolute conditions", "90-day shock"),
           modkind = c("Absolute", "Within-stratum"),
           panel = c("Minnesota", "Other flyway states", "Pooled"), sorted = FALSE)
res <- rbindlist(apply(grid, 1, function(g) {
  pn <- g[["panel"]]
  d <- if (pn == "Minnesota") cc[state == "Minnesota"]
       else if (pn == "Pooled") copy(cc) else cc[state != "Minnesota"]
  cat(sprintf("%-20s %-15s %-20s n=%d\n", g[["kind"]], g[["modkind"]], pn,
              sum(d$outbreak_binary))); flush.console()
  r <- try(run(d, pn, g[["modkind"]], g[["kind"]]), silent = TRUE)
  if (inherits(r, "try-error")) { cat("   failed:", attr(r, "condition")$message, "\n"); NULL } else r
}))
saveRDS(res, paste0(objects_folder, "si_survival_window.RDS"))
cat(sprintf("\n=== waterfowl effect (lag %d-%d d) by temperature over lag %d-%d d ===\n",
            min(BIRDWIN), max(BIRDWIN), min(SURV_WIN), max(SURV_WIN)))
cat(sprintf("    cold = %.1f degC, warm = %.1f degC. Prediction: cold > warm.\n\n", QS[1], QS[2]))
print(res, nrows = 100)
