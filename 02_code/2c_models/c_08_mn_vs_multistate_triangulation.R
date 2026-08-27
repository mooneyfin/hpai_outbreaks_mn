# c_08_mn_vs_multistate_triangulation.R
# Triangulation for the review response: does the Minnesota climate signal survive when
# (a) MN is analysed under the multi-state gridded design, and (b) the 6 other states are
# pooled in? Fits the identical 7-predictor DLNM to three daily subsets and two weekly
# subsets, and reports cumulative + early/late lag-window ORs (per +0.5 SD). If MN-only
# shows the climate signal but the pool does not -> heterogeneity (MN-specific). If MN-only
# is also null -> the original signal was a design/data artifact, not a robust effect.

rm(list = ls())
project.folder = paste0(print(here::here()), '/')
source(paste0(project.folder, 'create_folder_structure.R'))
source(paste0(functions.folder, 'script_initiate.R'))
conflict_prefer("year",  "lubridate", quiet = TRUE)
conflict_prefer("month", "lubridate", quiet = TRUE)

predictors <- c("runoff", "soil_moisture", "temperature", "precipitation",
                "wind_speed", "snow_cover", "anseriformes")
pred_labels <- c(runoff = "Runoff", soil_moisture = "Soil Moisture",
                 temperature = "Temperature", precipitation = "Precipitation",
                 wind_speed = "Wind Speed", snow_cover = "Snow Cover",
                 anseriformes = "Anseriformes")

# 1. Fit the primary DLNM (daily: met ns df3 + bird weekly-strata; weekly: ns df2).
fit_dlnm <- function(df, res) {
  setDT(df)
  if (res == "daily") { max_lag <- 28; arglag_met <- list(fun = "ns", df = 3)
    arglag_bird <- list(fun = "strata", breaks = c(7, 14, 21)) } else {
    max_lag <- 4;  arglag_met <- list(fun = "ns", df = 2); arglag_bird <- list(fun = "ns", df = 2) }
  cbs <- setNames(lapply(predictors, function(v) {
    X  <- as.matrix(df[, paste0(v, "_Lag", 0:max_lag), with = FALSE])
    al <- if (v == "anseriformes") arglag_bird else arglag_met
    crossbasis(X, lag = c(0, max_lag), argvar = list(fun = "lin"), arglag = al)
  }), paste0("cb_", predictors))
  list2env(cbs, environment())
  m <- clogit(as.formula(paste("outbreak_binary ~", paste(names(cbs), collapse = " + "),
                               "+ strata(stratum_id)")), method = "efron", data = df)
  list(model = m, cbs = cbs, max_lag = max_lag, n_cases = sum(df$outbreak_binary == 1))
}

# 2. OR (per +0.5 SD) over named lag windows, via coef+vcov contrasts.
window_or <- function(fit, windows) {
  beta <- coef(fit$model); V <- vcov(fit$model); lags <- 0:fit$max_lag
  rbindlist(lapply(predictors, function(v) {
    ix <- grep(paste0("^cb_", v, "v"), names(beta)); b <- beta[ix]; Vs <- V[ix, ix, drop = FALSE]
    lb <- do.call(dlnm::onebasis, c(list(x = lags), attr(fit$cbs[[paste0("cb_", v)]], "arglag")))
    rbindlist(lapply(names(windows), function(wn) {
      w <- windows[[wn]]; ct <- 0.5 * colSums(lb[lags >= w[1] & lags <= w[2], , drop = FALSE])
      lo <- sum(ct * b); se <- sqrt(as.numeric(t(ct) %*% Vs %*% ct))
      data.table(predictor = v, label = pred_labels[[v]], window = wn,
                 OR = exp(lo), lcl = exp(lo - 1.96 * se), ucl = exp(lo + 1.96 * se))
    }))
  }))
}

fmt <- function(or, lo, hi) sprintf("%.2f (%.2f,%.2f)%s", or, lo, hi, ifelse(lo > 1 | hi < 1, "*", " "))

# 3. Assemble subsets, fit, extract, and print side-by-side comparison per resolution.
o <- objects_folder
compare_res <- function(res, subsets, windows) {
  res_tab <- rbindlist(lapply(names(subsets), function(sn) {
    fit <- fit_dlnm(subsets[[sn]], res)
    wo  <- window_or(fit, windows)[, subset := sprintf("%s (n=%d)", sn, fit$n_cases)]
    wo
  }))
  res_tab[, cell := fmt(OR, lcl, ucl)]
  for (wn in names(windows)) {
    cat("\n#### ", toupper(res), " — ", wn, " lag window (OR per +0.5 SD; * = CI excludes 1)\n", sep = "")
    wide <- dcast(res_tab[window == wn], label ~ subset, value.var = "cell")
    setcolorder(wide, c("label", grep("=", names(wide), value = TRUE)))
    print(wide)
  }
  invisible(res_tab)
}

ms_daily  <- readRDS(paste0(o, "case_crossover_df_multistate_daily.RDS"))
ms_weekly <- readRDS(paste0(o, "case_crossover_df_multistate_weekly.RDS"))
setDT(ms_daily); setDT(ms_weekly)

daily_subsets <- list(
  `MN-original`   = readRDS(paste0(o, "case_crossover_df.RDS")),
  `MN-multistate` = ms_daily[state == "Minnesota"],
  `Pooled-7state` = ms_daily)
weekly_subsets <- list(
  `MN-multistate` = ms_weekly[state == "Minnesota"],
  `Pooled-7state` = ms_weekly)

daily_windows  <- list(cumulative = c(0, 28), `early(0-7d)` = c(0, 7), `late(21-28d)` = c(21, 28))
weekly_windows <- list(cumulative = c(0, 4),  `early(0-1wk)` = c(0, 1), `late(3-4wk)`  = c(3, 4))

cat("=================== DAILY ===================\n")
compare_res("daily",  daily_subsets,  daily_windows)
cat("\n=================== WEEKLY ==================\n")
compare_res("weekly", weekly_subsets, weekly_windows)
