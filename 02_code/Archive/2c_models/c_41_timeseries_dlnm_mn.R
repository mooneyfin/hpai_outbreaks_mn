############ MINNESOTA TIME-SERIES DLNM — 25-YEAR ANOMALY, ALL CELLS, NEGATIVE BINOMIAL ############
# A different DESIGN, not another parameterisation of the case-crossover.
#
# The case-crossover conditions on cells that had an outbreak and compares days within a cell, so
# it throws away every between-cell contrast by construction. This keeps all 3,334 Minnesota
# fishnet cells and models the daily count directly, so a cell-day where the climate was unusual
# is compared against cell-days where it was not, across space as well as time.
#
# The 25-year anomaly is the right exposure for this design. Absolute conditions would be almost
# entirely season, and a seasonal control term would then be fighting the exposure for the same
# variation. A departure from the 25-year weekly normal is close to season-free by construction.
#
# ABSOLUTE anomaly (x - normal), not z. Dividing by the weekly sd blows up for runoff and
# precipitation in cell-weeks that are historically almost always dry; see section 25 of the
# build plan.
#
# Negative binomial because the counts are overdispersed: aggregated over Minnesota the daily
# variance/mean is 2.12. Poisson is fitted alongside so the dispersion can be read off rather
# than assumed.

project.folder = paste0(print(here::here()), '/')
src <- readLines(paste0(project.folder, 'create_folder_structure.R'))
eval(parse(text = paste(src[!grepl("^\\s*rm\\(list", src)], collapse = "\n")))
source(paste0(functions.folder, 'script_initiate.R'))
suppressMessages({library(INLA); library(FNN); library(Matrix)})
source(paste0(functions.folder, 'inla_dlnm_helpers.R'))
INLA::inla.setOption(num.threads = 4)

NWEEK <- 4                                   # lag weeks 0-3, matching the weekly case-crossover
VARS  <- c("temperature", "precipitation", "soil_moisture", "wind_speed")
LAB   <- c(temperature = "Temperature", precipitation = "Precipitation",
           soil_moisture = "Soil moisture", wind_speed = "Wind speed",
           anseriformes = "Anseriformes")

# ---- 1. Minnesota fishnet panel ----
d <- setDT(readRDS(paste0(objects_folder, "fishnet_timeseries_daily_modeling.RDS")))
d[, date := as.Date(date)]
setorder(d, zone_id, date)
cat(sprintf("panel: %s cell-days, %d cells, %d days, %d outbreaks\n",
            format(nrow(d), big.mark = ","), uniqueN(d$zone_id), uniqueN(d$date),
            sum(d$outbreak_count)))

# ---- 2. attach the climatology by nearest centroid ----
# The fishnet and the ERA5 climatology are different grids with unrelated zone_id schemes
# (sequential integers vs "row,col" strings), so they have to be joined in space. The fishnet is
# 0.09 deg (~10 km) and the climatology is ERA5-Land native, so a nearest-centroid match is
# comfortably sub-cell: median 4.4 km, max 18.7 km.
MEAN <- paste0(VARS, "_mean")
clim <- rbindlist(lapply(
  list.files(env_data, pattern = "flyway8_climatology_1997_2021_w", full.names = TRUE),
  fread, select = c("zone_id", "state", "week_idx", MEAN),
  colClasses = list(character = "zone_id")))
clim <- clim[state == "Minnesota"][, state := NULL]

grid <- unique(d[, .(zone_id, lat, lon)])
hub  <- unique(clim[, .(zone_id)])
hubxy <- unique(rbindlist(lapply(
  list.files(env_data, pattern = "flyway8_climatology_1997_2021_w", full.names = TRUE),
  fread, select = c("zone_id", "state", "lat", "lon"),
  colClasses = list(character = "zone_id")))[state == "Minnesota", .(zone_id, lat, lon)])
kx <- cos(mean(grid$lat) * pi / 180)                 # crude equirectangular scaling, fine at this span
nn <- get.knnx(cbind(hubxy$lon * kx, hubxy$lat), cbind(grid$lon * kx, grid$lat), k = 1)
grid[, clim_id := hubxy$zone_id[nn$nn.index]]
cat(sprintf("climatology matched: median %.2f km, max %.2f km\n",
            median(nn$nn.dist) * 111, max(nn$nn.dist) * 111))

d <- merge(d, grid[, .(zone_id, clim_id)], by = "zone_id", all.x = TRUE)
d[, week_idx := pmin(as.integer(strftime(date, "%j")) %/% 7L, 51L)]
d <- merge(d, clim, by.x = c("clim_id", "week_idx"), by.y = c("zone_id", "week_idx"), all.x = TRUE)
setorder(d, zone_id, date)
for (v in VARS) set(d, j = paste0(v, "_an"), value = d[[v]] - d[[paste0(v, "_mean")]])
AN <- paste0(VARS, "_an")
cat(sprintf("anomalies complete for %.1f%% of cell-days\n", 100 * mean(!is.na(d$temperature_an))))

# ---- 3. weekly lag matrices ----
# Trailing 7-day means, shifted back a week at a time. shift() is within cell because the panel
# is ordered and grouped by zone_id, so lags never bleed across cell boundaries.
roll7 <- function(x) data.table::frollmean(x, 7, align = "right")
for (v in c(AN, "anseriformes")) {
  d[, (paste0(v, "_w0")) := roll7(get(v)), by = zone_id]
  for (w in 1:(NWEEK - 1))
    d[, (paste0(v, "_w", w)) := shift(get(paste0(v, "_w0")), 7L * w), by = zone_id]
}
need <- as.vector(outer(c(AN, "anseriformes"), paste0("_w", 0:(NWEEK - 1)), paste0))
d <- d[complete.cases(d[, ..need])]
cat(sprintf("after lag construction: %s cell-days, %d outbreaks\n",
            format(nrow(d), big.mark = ","), sum(d$outbreak_count)))

# ---- 4. crossbasis per exposure ----
# Same shapes as the weekly case-crossover so the two designs are comparable: ns(2) on the
# exposure for meteorology, linear for waterfowl, ns(2) across the lag axis.
cb_of <- function(v, lin = FALSE) {
  m <- as.matrix(d[, paste0(v, "_w", 0:(NWEEK - 1)), with = FALSE])
  m <- scale(m)                                   # one centre/scale for the whole lag matrix
  cb <- dlnm::crossbasis(m, lag = c(0, NWEEK - 1),
                         argvar = if (lin) list(fun = "lin") else list(fun = "ns", df = 2),
                         arglag = list(fun = "ns", df = 2))
  colnames(cb) <- paste0("cb_", v, ".", colnames(cb))
  cb
}
bases <- c(lapply(setNames(AN, AN), cb_of), list(anseriformes = cb_of("anseriformes", lin = TRUE)))
# unname() matters: cbind on a NAMED list of data frames prefixes every column with the list
# name, so cb_temperature_an.v1.l1 would silently become temperature_an.cb_temperature_an.v1.l1
X <- do.call(cbind, unname(lapply(bases, function(b) as.data.frame(unclass(b)))))

# ---- 5. model frame ----
# cell: iid, soaks up everything time-invariant about a cell that we have not measured
# t   : rw2 on sequential day, the temporal control
# log1p(feedlot_count) as a covariate rather than an offset: 1,189 cells have no registered
#       feedlot at all, and an offset would either drop them or assert they carry one farm's
#       worth of risk. As a covariate the data sets the slope and every cell stays in.
md <- data.frame(y = d$outbreak_count, X, check.names = FALSE,
                 cell = as.integer(as.factor(d$zone_id)),
                 t    = as.integer(d$date - min(d$date)) + 1L,
                 doy  = as.integer(d$doy),
                 lfeed = log1p(d$feedlot_count))
cbn <- unlist(lapply(bases, colnames), use.names = FALSE)
stopifnot(all(cbn %in% names(md)))
PC  <- list(prec = list(prior = "pc.prec", param = c(1, 0.01)))

fit_ts <- function(tcontrol, family = "nbinomial") {
  f <- as.formula(paste("y ~", paste(cbn, collapse = " + "), "+ lfeed +",
                        "f(cell, model = 'iid', hyper = PC) +", tcontrol))
  INLA::inla(f, family = family, data = md,
             control.fixed = list(prec = 4, prec.intercept = 0.01),
             control.compute = list(waic = TRUE, config = TRUE),
             control.predictor = list(compute = FALSE),
             control.inla = list(strategy = "gaussian", int.strategy = "eb"))
}

# rate ratio per +0.5 SD summed over a lag window, read off the fitted crossbasis exactly as
# the case-crossover does it
rr <- function(fit, v, w) {
  cb <- bases[[v]]
  bs <- function(x, a) do.call(dlnm::onebasis, c(list(x = x), a))
  ev <- as.numeric(bs(0.5, attr(cb, "argvar")) - bs(0, attr(cb, "argvar")))
  lv <- colSums(bs(w, attr(cb, "arglag")), na.rm = TRUE)
  wt <- as.vector(t(outer(ev, lv)))
  cv <- dlnm_coef_vcov(fit, colnames(cb))       # joint covariance, off config when needed
  lr <- sum(wt * cv$coef); se <- sqrt(as.numeric(t(wt) %*% cv$vcov %*% wt))
  sprintf("%.2f (%.2f, %.2f)", exp(lr), exp(lr - 1.96 * se), exp(lr + 1.96 * se))
}

MODELS <- list(
  `NB, RW2 on day`       = list(t = "f(t, model = 'rw2', hyper = PC)",     fam = "nbinomial"),
  `NB, seasonal spline`  = list(t = "f(doy, model = 'rw2', cyclic = TRUE, hyper = PC)", fam = "nbinomial"),
  `Poisson, RW2 on day`  = list(t = "f(t, model = 'rw2', hyper = PC)",     fam = "poisson"))

res <- list()
for (nm in names(MODELS)) {
  cat(sprintf("\n--- %s ---\n", nm))
  fit <- try(fit_ts(MODELS[[nm]]$t, MODELS[[nm]]$fam), silent = TRUE)
  if (inherits(fit, "try-error")) { cat("failed:", attr(fit, "condition")$message, "\n"); next }
  cat(sprintf("WAIC %.1f\n", fit$waic$waic))
  if (MODELS[[nm]]$fam == "nbinomial") {
    h <- fit$summary.hyperpar[grep("size", rownames(fit$summary.hyperpar)), ]
    if (nrow(h)) cat(sprintf("NB size (dispersion): %.2f (%.2f, %.2f)\n",
                             h[1, "mean"], h[1, "0.025quant"], h[1, "0.975quant"]))
  }
  res[[nm]] <- rbindlist(lapply(names(bases), function(v) rbindlist(lapply(0:(NWEEK - 1),
    function(w) data.table(Model = nm, Exposure = unname(LAB[sub("_an$", "", v)]),
                           `Lag week` = paste0("Week ", w), est = rr(fit, v, w))))))
  res[[nm]][, WAIC := round(fit$waic$waic, 1)]
}
out <- rbindlist(res)
saveRDS(out, paste0(objects_folder, "si_timeseries_dlnm_mn.RDS"))
cat("\n=== Minnesota time-series DLNM, 25-year absolute anomaly ===\n")
print(dcast(out, Exposure + Model ~ `Lag week`, value.var = "est"), nrows = 100)
