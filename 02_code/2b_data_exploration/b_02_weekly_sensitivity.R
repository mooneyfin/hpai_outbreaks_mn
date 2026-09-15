############ SUPPLEMENTARY — WEEKLY EXPOSURE RESOLUTION ############
# The original submission modelled exposures at weekly resolution; this revision moved to
# daily. This is the bridge: the same design, the same case definition, the same strata, but
# exposures binned into weekly means over a 4-week lag rather than 29 daily lags.
#
# Lag bins are days 0-6, 7-13, 14-20, 21-27, so the window matches the daily model's 0-28 d
# and a reader can compare like with like. Binning drops the crossbasis from 29 lag columns
# to 4, which is most of the reason the weekly model was ever attractive at this sample size.

rm(list = ls())
project.folder = paste0(print(here::here()), '/')
source(paste0(project.folder, 'create_folder_structure.R'))
source(paste0(functions.folder, 'script_initiate.R'))
suppressMessages({library(INLA); library(Matrix); library(dlnm); library(flextable); library(officer)})
INLA::inla.setOption(num.threads = 1)
source(paste0(functions.folder, 'inla_dlnm_helpers.R'))

PREC  <- 4
NWEEK <- 4                       # lag weeks 0-3, i.e. days 0-27
SH    <- c("temperature_shock", "soil_moisture_shock", "runoff_shock", "wind_speed_shock")
LV    <- c("temperature", "soil_moisture", "runoff", "wind_speed")
BIRD  <- "anseriformes"
LABS  <- c(temperature_shock = "Temperature", soil_moisture_shock = "Soil moisture",
           runoff_shock = "Runoff", wind_speed_shock = "Wind speed",
           temperature = "Temperature", soil_moisture = "Soil moisture",
           runoff = "Runoff", wind_speed = "Wind speed", anseriformes = "Anseriformes")

D <- setDT(readRDS(paste0(objects_folder, "case_crossover_df_shock_post7.RDS")))
if (anyDuplicated(names(D))) D <- D[, which(!duplicated(names(D))), with = FALSE]

# 1. Weekly bins from the daily lag columns, then z-scored the same way the daily ones were
week_mat <- function(dat, v) {
  m <- vapply(0:(NWEEK - 1), function(w) {
    cols <- paste0(v, "_Lag", (w * 7):(w * 7 + 6))
    rowMeans(as.matrix(dat[, ..cols]), na.rm = TRUE)
  }, numeric(nrow(dat)))
  apply(m, 2, function(x) as.numeric(scale(x)))
}

# 2. Fit by hand rather than through assemble_cc_dlnm, because the lag axis is weeks here
fit_weekly <- function(dat, vars) {
  bases <- list(); cb_cols <- list()
  dd <- data.frame(outbreak_binary = dat$outbreak_binary)
  for (v in vars) {
    av <- if (v == BIRD) list(fun = "lin") else list(fun = "ns", df = 2)
    cb <- dlnm::crossbasis(week_mat(dat, v), lag = c(0, NWEEK - 1),
                           argvar = av, arglag = list(fun = "ns", df = 2))
    colnames(cb) <- paste0("cb_", v, ".", colnames(cb))
    bases[[v]] <- cb; cb_cols[[v]] <- colnames(cb)
    dd <- cbind(dd, as.data.frame(unclass(cb)))
  }
  dd$id_stratum <- as.integer(as.factor(dat$stratum_id))
  o <- list(data = dd, bases = bases, cb_cols = cb_cols, met_vars = vars,
            contrast_steps = list(), weekly = list(), weekly_vars = character(0),
            max_lag = NWEEK - 1, max_lag_by_var = list(),
            argvar = list(fun = "ns", df = 2), argvar_by_var = list(),
            arglag = list(fun = "ns", df = 2), arglag_by_var = list())
  list(obj = o, fit = suppressWarnings(fit_cc_inla_dlnm(o, fixed_prec = PREC)))
}

cell <- function(e) sprintf("%.2f (%.2f, %.2f)", e[["OR"]], e[["low"]], e[["high"]])
res <- rbindlist(lapply(c("Minnesota", "Pooled"), function(pn) {
  d <- if (pn == "Minnesota") D[state == "Minnesota"] else copy(D)
  rbindlist(lapply(c(anomaly = "anomaly", absolute = "absolute"), function(kind) {
    vars <- c(if (kind == "anomaly") SH else LV, BIRD)
    r <- fit_weekly(d, vars)
    cat(sprintf("%-11s %-8s n = %3d | WAIC %.1f\n", pn, kind, sum(d$outbreak_binary),
                r$fit$waic$waic))
    rbindlist(lapply(vars, function(v) {
      wk <- vapply(0:(NWEEK - 1), function(w)
        cell(met_lag_effect(r$fit, r$obj, v, w)), character(1))
      cu <- cell(met_lag_effect(r$fit, r$obj, v, 0:(NWEEK - 1)))
      as.data.table(setNames(as.list(c(pn, kind, unname(LABS[v]), wk, cu)),
        c("Panel", "Parameterisation", "Predictor",
          paste0("Lag week ", 0:(NWEEK - 1)), "Cumulative 0-4 weeks")))
    }))
  }))
}))

saveRDS(res, paste0(objects_folder, "tableS6_weekly_sensitivity.RDS"))
is_sig <- function(x) {
  m <- regmatches(x, regexec("\\(([0-9.]+), ([0-9.]+)\\)", x))
  vapply(m, function(z) length(z) == 3 && (as.numeric(z[2]) > 1 | as.numeric(z[3]) < 1), TRUE)
}
ft <- flextable(res) |>
  add_header_row(values = c("", "Odds ratio (95% credible interval) per +0.5 SD"),
                 colwidths = c(3, ncol(res) - 3)) |>
  fontsize(size = 10, part = "all") |> autofit() |>
  align(align = "center", part = "all") |> theme_booktabs()
for (j in names(res)[-(1:3)]) {
  hit <- which(is_sig(res[[j]])); if (length(hit)) ft <- bold(ft, i = hit, j = j)
}
read_docx() |> body_add_flextable(ft) |>
  print(target = file.path(tables_main_folder, "TableS6_weekly_sensitivity.docx"))

cat("\n=== Weekly exposure resolution, 4-week lag ===\n")
print(res[Parameterisation == "anomaly", -"Parameterisation"])
cat("\nwrote TableS6_weekly_sensitivity.docx\n")
