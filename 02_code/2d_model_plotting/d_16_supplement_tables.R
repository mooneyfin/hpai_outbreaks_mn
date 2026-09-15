############ SUPPLEMENT — S2 TO S6 ############
# S1 is the colleague's WQS mixture and isn't built here.
#
#   S2  Primary results in full - both models, both estimate types, all three panels
#   S3  Referent design
#   S4  Exposure set, including the three runoff tests
#   S5  Lag structure and prior
#   S6  Sensitivity: climatological anomaly, weekly bins, cluster-robust intervals
#   S7  Season interaction
#
# S2 used to be Table 2 in the main text. It moved because Figure 3 already shows every
# cumulative estimate it carried - same five exposures, same four windows, same three panels,
# both models - so the table was the figure again in numbers. What it uniquely adds is the
# lag-specific column, and per-day contributions sitting around 1.03 / 0.97 are a lookup, not a
# main-text exhibit. The three significant estimates go in the Results prose instead, which is
# where anyone extracting numbers will actually look.
#
# TWO PRIMARY MODELS, and every table carries both: absolute conditions and the 90-day rolling
# anomaly. They disagree about soil moisture (1.05 vs 0.81) and that disagreement is a result,
# not a footnote - it's the clearest evidence we have about which associations survive taking
# the season out.
#
# Two things the tables must not imply:
#
# 1. A dWAIC COLUMN IS NOT VALID EVERYWHERE. Changing the referent design changes which
#    observations exist, so those WAICs sit on different data and differencing them just rewards
#    whichever design deletes the most rows. Those cells are empty on purpose. Same for the
#    prior, where WAIC falls monotonically as you shrink toward the null - that's shrinkage, not
#    fit, so the prior is chosen on its implied OR range instead.
#
# 2. The "% of annual range" number on the design rows IS NOT A NAMED OR VALIDATED DIAGNOSTIC.
#    It's arithmetic on the 25-year climatology - how far the climatological mean temperature
#    moves across the days inside one stratum. Reported in °C; these are temperature SPREADS, so the degC and K values are numerically identical, but the paper reports temperature in Celsius throughout. Report it as description, never as a criterion.
#    The design argument rests on the literature: Janes, Sheppard & Lumley (2005) for why a
#    fixed calendar partition keeps the conditional likelihood unbiased, and Bateson & Schwartz
#    (1999) for narrower referent windows controlling season better.

rm(list = ls())
project.folder = paste0(print(here::here()), '/')
source(paste0(project.folder, 'create_folder_structure.R'))
source(paste0(functions.folder, 'script_initiate.R'))
suppressMessages({library(INLA); library(Matrix); library(flextable); library(officer)})
INLA::inla.setOption(num.threads = 1)
source(paste0(functions.folder, 'inla_dlnm_helpers.R'))

PREC <- 4; NS2 <- list(fun = "ns", df = 2); LIN <- list(fun = "lin")
SH   <- c("temperature_shock", "precipitation_shock", "soil_moisture_shock", "wind_speed_shock")
BIRD <- "anseriformes"
cell <- function(e) sprintf("%.2f (%.2f, %.2f)", e[["OR"]], e[["low"]], e[["high"]])
blank <- function(n) rep("", n)

D <- setDT(readRDS(paste0(objects_folder, "case_crossover_df_shock_post7.RDS")))
if (anyDuplicated(names(D))) D <- D[, which(!duplicated(names(D))), with = FALSE]

# one fit on the pooled panel, reporting EVERY exposure. the table has to let a reader see the
# whole result set move (or not) under each alternative; showing two terms invites the obvious
# question of why those two.
ALL  <- c(SH, BIRD)
WIN  <- list(temperature_shock = 0:7, precipitation_shock = 0:7, soil_moisture_shock = 0:7,
             wind_speed_shock = 0:7, runoff_shock = 0:7, anseriformes = 22:28,
             temperature = 0:7, precipitation = 0:7, soil_moisture = 0:7,
             wind_speed = 0:7, runoff = 0:7)
LAB  <- c(temperature_shock = "Temperature 0-7 d", precipitation_shock = "Precipitation 0-7 d",
          soil_moisture_shock = "Soil moisture 0-7 d", wind_speed_shock = "Wind speed 0-7 d",
          runoff_shock = "Runoff 0-7 d", anseriformes = "Anseriformes 22-28 d",
          # absolute-scale twins, same column labels so the design block lines up with the rest
          temperature = "Temperature 0-7 d", precipitation = "Precipitation 0-7 d",
          soil_moisture = "Soil moisture 0-7 d", wind_speed = "Wind speed 0-7 d",
          runoff = "Runoff 0-7 d")
LEVEL <- c("temperature", "precipitation", "soil_moisture", "wind_speed")
COLS <- unname(LAB[c(SH, "runoff_shock", BIRD)])

est_row <- function(f, o, mv, max_lag = 28) {
  out <- setNames(as.list(rep("--", length(COLS))), COLS)     # "--" = not in this model
  for (v in names(LAB)) {
    if (!v %in% mv) next
    w <- WIN[[v]]; w <- w[w <= max_lag]
    if (!length(w)) next
    out[[LAB[[v]]]] <- cell(met_lag_effect(f, o, v, w))
  }
  out
}

fit_spec <- function(vars, bird = TRUE, max_lag = 28, lag_df = 2, dat = D) {
  mv <- if (bird) c(vars, BIRD) else vars
  av <- c(setNames(rep(list(NS2), length(vars)), vars),
          if (bird) setNames(list(LIN), BIRD))
  o <- assemble_cc_dlnm(dat, met_vars = mv, weekly_vars = character(0), max_lag = max_lag,
                        argvar_by_var = av, arglag = list(fun = "ns", df = lag_df))
  f <- suppressWarnings(fit_cc_inla_dlnm(o, fixed_prec = PREC))
  c(list(waic = f$waic$waic), est_row(f, o, mv, max_lag))
}

mkrow <- function(domain, choice, r, diag = "", refs = "", dwaic = "") {
  cbind(data.table(Domain = domain, Choice = choice, `Referents per case` = refs,
                   Diagnostic = diag, dWAIC = dwaic),
        as.data.table(r[COLS]))
}


# ---- shared: the two primary models ----
MODELS <- list(`Absolute` = LEVEL, `90-day anomaly` = SH)
RUNOFF <- c(Absolute = "runoff", `90-day anomaly` = "runoff_shock")

refw <- vapply(names(MODELS), function(m) fit_spec(MODELS[[m]], dat = D)$waic, numeric(1))
dw   <- function(w, m) sprintf("%+.1f", w - refw[[m]])
cat(sprintf("reference WAIC | absolute %.1f | anomaly %.1f\n", refw[[1]], refw[[2]]))

# every table is built the same way: loop the two models, loop that table's variants.
build <- function(variants, fit_one, extra = function(nm, m, r) list()) {
  rbindlist(lapply(names(MODELS), function(m) {
    rbindlist(lapply(names(variants), function(nm) {
      r <- fit_one(variants[[nm]], m)
      cbind(data.table(Model = m, Choice = nm), as.data.table(extra(nm, m, r)),
            as.data.table(r[COLS]))
    }))
  }))
}

# ---- S3 Referent design ----
# No dWAIC: each row is a different set of observations. The diagnostic column is plain
# description (see header) - how much of the annual temperature cycle survives inside a stratum.
ex <- setDT(readRDS(paste0(objects_folder, "daily_shock_predictors.RDS")))
ex[, date := as.Date(date)]
SHOCK <- setdiff(names(ex), c("zone_id", "date"))
attach_shocks <- function(f) {
  cc <- setDT(readRDS(paste0(objects_folder, f)))
  if (anyDuplicated(names(cc))) cc <- cc[, which(!duplicated(names(cc))), with = FALSE]
  cc[, date := as.Date(date)]
  for (L in 0:28) {
    sr <- ex[, c(.(zone_id = zone_id, date = date + L), .SD), .SDcols = SHOCK]
    setnames(sr, SHOCK, paste0(SHOCK, "_Lag", L))
    cc <- merge(cc, sr, by = c("zone_id", "date"), all.x = TRUE, sort = FALSE)
  }
  for (v in SHOCK) for (L in 0:28) {
    cn <- paste0(v, "_Lag", L); set(cc, j = cn, value = as.numeric(scale(cc[[cn]])))
  }
  need <- paste0(rep(SHOCK, each = 29), "_Lag", 0:28)
  cc <- cc[complete.cases(cc[, ..need])]
  cc[, `:=`(nc = sum(outbreak_binary), nr = sum(outbreak_binary == 0)), by = stratum_id]
  cc[nc == 1L & nr > 0L]
}

DESIGNS <- list(
  `Symmetric bidirectional, +/-14 to 28 d` = list(f = "case_crossover_df_sym_g14.RDS",              s = "55.4 d, 4.74 °C (45%)"),
  `Time-stratified, two-month strata`      = list(f = "case_crossover_df_timestrat_bimonth_nowo.RDS", s = "59.1 d, 3.59 °C (34%)"),
  `Time-stratified, one-month strata`      = list(f = "case_crossover_df_timestrat_month.RDS",      s = "29.4 d, 2.00 °C (19%)"),
  `One-month strata, 7-day post exclusion  [selected]` = list(f = NA,                                       s = "28.6 d, 1.98 °C (19%)"))
DPANEL <- new.env()
get_panel <- function(f) {
  if (is.na(f)) return(D)
  if (is.null(DPANEL[[f]])) assign(f, attach_shocks(f), envir = DPANEL)
  DPANEL[[f]]
}

S3 <- build(DESIGNS,
  function(g, m) fit_spec(MODELS[[m]], dat = get_panel(g$f)),
  function(nm, m, r) {
    d <- get_panel(DESIGNS[[nm]]$f)
    list(`Referents per case` = sprintf("%.1f", sum(d$outbreak_binary == 0) / sum(d$outbreak_binary)),
         `Season left in a stratum` = DESIGNS[[nm]]$s)
  })

# ---- S2 Primary results in full ----
# Straight off the cached c_17 fits. Cumulative AND lag-specific, because a DLNM paper should
# show the per-day contribution as well as what it accumulates to - it just doesn't need to do
# that in the main text.
pf   <- readRDS(paste0(objects_folder, "primary_shock_level_fits.RDS"))
KEY  <- c(Minnesota = "Minnesota", `Other flyway states` = "Other flyway", Pooled = "Pooled")
KIND <- c(`Absolute` = "level", `90-day anomaly` = "shock")
EST  <- c(Cumulative = "cumulative", `Lag-specific` = "lagspecific")

S2 <- rbindlist(lapply(names(KIND), function(mn)
  rbindlist(lapply(names(KEY), function(pn) {
    r <- pf[[paste(KEY[[pn]], KIND[[mn]], sep = "|")]]
    rbindlist(lapply(names(EST), function(en)
      cbind(data.table(Model = mn, Panel = pn, Estimate = en), copy(r[[EST[[en]]]]))),
      fill = TRUE)
  }))), fill = TRUE)

# the two estimate types don't share columns - cumulative reports four windows (0-7 ... 22-28),
# lag-specific reports five anchor days (0, 7, 14, 21, 28) - so every row is blank in the half
# that doesn't apply to it. fill = TRUE leaves those as NA, which prints as "NA" in the docx.
WCOL <- names(INLA_LAG_WINDOWS)
ACOL <- paste0("Lag ", INLA_LAG_ANCHORS, " days")
setcolorder(S2, c("Model", "Panel", "Estimate", "Predictor", WCOL, ACOL))
for (j in c(WCOL, ACOL)) set(S2, which(is.na(S2[[j]])), j, "")

# ---- S4 Exposure set, with runoff tested three ways ----
# Runoff was the strongest predictor in the original submission, so the supplement has to show
# what it reads under the current spec rather than just assert it was dropped.
exp_variants <- function(m) {
  b <- MODELS[[m]]; ro <- RUNOFF[[m]]           # b = temp, precip, soil, wind in that order
  list(`Temperature + precipitation + soil moisture + wind + waterfowl  [selected]` = list(v = b, bird = TRUE),
       `drop temperature`   = list(v = b[-1], bird = TRUE),
       `drop precipitation` = list(v = b[-2], bird = TRUE),
       `drop soil moisture` = list(v = b[-3], bird = TRUE),
       `drop wind speed`    = list(v = b[-4], bird = TRUE),
       `drop waterfowl`     = list(v = b,     bird = FALSE),
       `add runoff to the selected model` = list(v = c(b, ro),      bird = TRUE),
       `runoff in place of precipitation` = list(v = c(b[-2], ro),  bird = TRUE),
       `runoff as the only water term`    = list(v = c(b[c(1, 4)], ro), bird = TRUE))
}
S4 <- rbindlist(lapply(names(MODELS), function(m) {
  vs <- exp_variants(m)
  rbindlist(lapply(names(vs), function(nm) {
    r <- fit_spec(vs[[nm]]$v, bird = vs[[nm]]$bird)
    cat(sprintf("  S4 %-14s | %-62s dWAIC %6s\n", m, nm, dw(r$waic, m)))
    cbind(data.table(Model = m, Choice = nm, dWAIC = dw(r$waic, m)), as.data.table(r[COLS]))
  }))
}))

# ---- S5 Lag structure and prior ----
# Lag window stops at 28 d: 42 and 56 need lag columns this panel doesn't carry, so they'd sit
# on fewer rows and stop being differenceable against the rest.
LAGS <- list(`Lag 0-28 d, ns(df = 2)  [selected]` = list(m = 28, d = 2),
             `Lag 0-14 d, ns(df = 2)` = list(m = 14, d = 2),
             `Lag 0-21 d, ns(df = 2)` = list(m = 21, d = 2),
             `Lag 0-28 d, ns(df = 3)` = list(m = 28, d = 3),
             `Lag 0-28 d, ns(df = 4)` = list(m = 28, d = 4))
s5a <- rbindlist(lapply(names(MODELS), function(m) {
  rbindlist(lapply(names(LAGS), function(nm) {
    r <- fit_spec(MODELS[[m]], max_lag = LAGS[[nm]]$m, lag_df = LAGS[[nm]]$d)
    cat(sprintf("  S5 %-14s | %-62s dWAIC %6s\n", m, nm, dw(r$waic, m)))
    cbind(data.table(Model = m, Choice = nm, dWAIC = dw(r$waic, m), Diagnostic = ""),
          as.data.table(r[COLS]))
  }))
}))

PRIORS <- c(0.001, 1, 4, 16, 64, 256)
set.seed(1)
s5b <- rbindlist(lapply(names(MODELS), function(m) {
  mv <- c(MODELS[[m]], BIRD)
  o <- assemble_cc_dlnm(D, met_vars = mv, weekly_vars = character(0), max_lag = 28,
                        argvar_by_var = c(setNames(rep(list(NS2), length(MODELS[[m]])), MODELS[[m]]),
                                          setNames(list(LIN), BIRD)),
                        arglag = list(fun = "ns", df = 2))
  cb <- o$bases[[MODELS[[m]][3]]]                      # soil moisture, either scale
  ob <- function(x, a) do.call(dlnm::onebasis, c(list(x = x), a))
  vv <- as.numeric(ob(0.5, attr(cb, "argvar")) - ob(0, attr(cb, "argvar")))
  cvec <- as.vector(t(outer(vv, colSums(ob(0:7, attr(cb, "arglag")), na.rm = TRUE))))
  rbindlist(lapply(PRIORS, function(pr) {
    f  <- suppressWarnings(fit_cc_inla_dlnm(o, fixed_prec = pr))
    OR <- exp(as.numeric(matrix(rnorm(20000 * length(cvec), 0, 1 / sqrt(pr)),
                                ncol = length(cvec)) %*% cvec))
    lab <- sprintf("Prior N(0, %.3g^2)%s", 1 / sqrt(pr), if (pr == 4) "  [selected]" else "")
    r <- c(list(waic = f$waic$waic), est_row(f, o, mv))
    cat(sprintf("  S5 %-14s | %-62s dWAIC %6s\n", m, lab, dw(r$waic, m)))
    cbind(data.table(Model = m, Choice = lab, dWAIC = dw(r$waic, m),
                     Diagnostic = sprintf("implied OR %.2f-%.2f",
                                          quantile(OR, .025), quantile(OR, .975))),
          as.data.table(r[COLS]))
  }))
}))
S5 <- rbind(s5a, s5b)

# ---- S6 Sensitivity: what's left over ----
tp <- setDT(readRDS(paste0(objects_folder, "si_three_parameterisations.RDS")))
s6a <- tp[panel == "Pooled" & Predictor %in% c("Soil moisture", "Temperature", "Anseriformes"),
          .(Block = "", Variant = Parameterisation, Predictor,
            `Lag 0-7 d` = `Lag 0-7 days`, `Lag 22-28 d` = `Lag 22-28 days`,
            Note = paste0("n = ", n))]
s6a[1, Block := "Exposure parameterisation"]

wk <- setDT(readRDS(paste0(objects_folder, "tableS6_weekly_sensitivity.RDS")))
s6b <- wk[Panel == "Pooled" & Predictor %in% c("Soil moisture", "Temperature", "Anseriformes"),
          .(Block = "", Variant = Parameterisation, Predictor,
            `Lag 0-7 d` = `Lag week 0`, `Lag 22-28 d` = `Lag week 3`,
            Note = "weekly bins, 4-week lag")]
s6b[1, Block := "Temporal resolution"]

cl <- setDT(readRDS(paste0(objects_folder, "si_followup_sensitivities.RDS"))$cluster)
s6c <- cl[panel == "Pooled", .(Block = "", Variant = paste("Cluster bootstrap on cells,", window),
                               Predictor, `Lag 0-7 d` = model_based, `Lag 22-28 d` = cluster_robust,
                               Note = sprintf("interval width x%.2f", width_ratio))]
s6c[1, Block := "Clustered intervals"]
S6 <- rbind(s6a, s6b, s6c)

# ---- S7 Season interaction ----
S7 <- setDT(readRDS(paste0(objects_folder, "si_season_interaction.RDS")))[
  panel == "Pooled", .(Model = parameterisation, Exposure = exposure, Window = window,
                       Spring = spring, `Outside spring` = outside, Ratio = ratio,
                       `p interaction` = p_interaction)]

TABLES <- list(S2 = S2, S3 = S3, S4 = S4, S5 = S5, S6 = S6, S7 = S7)
saveRDS(TABLES, paste0(objects_folder, "supplement_tables.RDS"))

is_sig <- function(x) {
  m <- regmatches(x, regexec("\\(([0-9.]+), ([0-9.]+)\\)", x))
  vapply(m, function(z) length(z) == 3 && (as.numeric(z[2]) > 1 | as.numeric(z[3]) < 1), TRUE)
}
for (nm in names(TABLES)) {
  tt <- TABLES[[nm]]
  ft <- flextable(tt) |> fontsize(size = 8, part = "all") |> autofit() |>
    align(align = "center", part = "all") |> theme_booktabs()
  for (j in names(tt)) {
    if (!is.character(tt[[j]]) || !any(grepl("^[0-9.]+ \\(", tt[[j]]))) next
    hit <- which(is_sig(tt[[j]])); if (length(hit)) ft <- bold(ft, i = hit, j = j)
  }
  read_docx() |> body_add_flextable(ft) |>
    print(target = file.path(tables_main_folder, sprintf("Table%s.docx", nm)))
  cat(sprintf("\n===== %s (%d rows) =====\n", nm, nrow(tt))); print(tt, nrows = 60)
}
cat("\nwrote TableS2-S7 (.docx + supplement_tables.RDS)\n")
