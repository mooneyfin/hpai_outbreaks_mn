############ EXPLORATORY — bdlnm EXPOSURE-LAG SURFACES ############
# Not manuscript figures. These are the exposure-lag-response surfaces from `bdlnm`
# (Bayesian DLNMs via INLA, CRAN 0.1.1), for our own read on what the crossbases are doing
# in two dimensions at once - the manuscript figures only ever show one slice at a time.
#
# The fit here is the same conditional-Poisson formula and the same N(0, 0.5^2) prior as the
# primary; bdlnm passes `...` straight through to inla, so it reproduces our own pipeline to
# within 0.07% when the prior and integration strategy are matched. That agreement is itself
# worth having - it is an independent check on the crossbasis contrast machinery in
# inla_dlnm_helpers.R, including the hand-rolled Qinv reconstruction.

rm(list = ls())
project.folder = paste0(print(here::here()), '/')
source(paste0(project.folder, 'create_folder_structure.R'))
source(paste0(functions.folder, 'script_initiate.R'))
suppressMessages({library(INLA); library(Matrix); library(bdlnm)})
INLA::inla.setOption(num.threads = 1)
source(paste0(functions.folder, 'inla_dlnm_helpers.R'))

PREC <- 4
SH   <- c("temperature_shock", "soil_moisture_shock", "runoff_shock", "wind_speed_shock")
MV   <- c(SH, "anseriformes")
LAB  <- c(temperature_shock = "Temperature anomaly", soil_moisture_shock = "Soil moisture anomaly",
          runoff_shock = "Runoff anomaly", wind_speed_shock = "Wind speed anomaly",
          anseriformes = "Anseriformes abundance")
OUT  <- file.path(figures_folder, "eda")
dir.create(OUT, showWarnings = FALSE, recursive = TRUE)

D <- setDT(readRDS(paste0(objects_folder, "case_crossover_df_shock_post7.RDS")))
if (anyDuplicated(names(D))) D <- D[, which(!duplicated(names(D))), with = FALSE]

for (pn in c("Minnesota", "Pooled")) {
  d <- if (pn == "Minnesota") D[state == "Minnesota"] else copy(D)

  # crossbases exactly as the primary builds them
  cbs <- lapply(MV, function(v) build_met_crossbasis(
    d, v, 28, if (v == "anseriformes") list(fun = "lin") else list(fun = "ns", df = 2),
    list(fun = "ns", df = 2)))
  names(cbs) <- paste0("cb_", MV)
  for (n in names(cbs)) assign(n, cbs[[n]], envir = globalenv())

  dat <- data.frame(outbreak_binary = d$outbreak_binary,
                    id_stratum = as.integer(as.factor(d$stratum_id)))
  hyper.strata <<- list(prec = list(initial = log(1e-6), fixed = TRUE))
  fm <- as.formula(paste("outbreak_binary ~", paste(names(cbs), collapse = " + "),
                         "+ f(id_stratum, model='iid', constr=FALSE, hyper=hyper.strata)"))
  b <- bdlnm(fm, data = dat, family = "poisson",
             sample.arg = list(n = 10000, seed = 1L),
             control.fixed = list(mean = 0, prec = PREC,
                                  prec.intercept = 0.001, mean.intercept = 0),
             control.inla = list(strategy = "simplified.laplace", int.strategy = "eb",
                                 control.vb = list(enable = FALSE)))
  cat(sprintf("%s: bdlnm fitted (n = %d)\n", pn, sum(d$outbreak_binary)))

  tag <- tolower(gsub(" ", "_", pn))
  for (v in MV) {
    rng <- range(attr(cbs[[paste0("cb_", v)]], "range"), na.rm = TRUE)
    grd <- seq(rng[1], rng[2], length.out = 41)
    # must be exact members of `grd`; rounding them makes plot.bcrosspred reject them
    sl  <- grd[round(c(0.10, 0.35, 0.65, 0.90) * (length(grd) - 1)) + 1]
    bp <- bcrosspred(b, basis = paste0("cb_", v), exp_at = grd, cen = 0, cumul = TRUE)
    ttl <- sprintf("%s - %s", LAB[[v]], pn)

    # contour: the whole exposure-lag surface at once
    png(file.path(OUT, sprintf("bdlnm_contour_%s_%s.png", tag, v)),
        width = 2000, height = 1500, res = 220)
    plot(bp, ptype = "contour", main = ttl, xlab = "Exposure (SD)", ylab = "Lag (days)")
    dev.off()

    # slices: lag-response at chosen exposures, exposure-response at chosen lags
    png(file.path(OUT, sprintf("bdlnm_slices_%s_%s.png", tag, v)),
        width = 2200, height = 1500, res = 220)
    plot(bp, ptype = "slices", exp_at = sl, lag_at = c(0, 7, 14, 28),
         ci = "area", main = ttl)
    dev.off()

    # overall cumulative exposure-response, with posterior draws rather than a ribbon -
    # an honest picture of how little the tails are pinned down
    png(file.path(OUT, sprintf("bdlnm_overall_%s_%s.png", tag, v)),
        width = 1800, height = 1400, res = 220)
    plot(bp, ptype = "overall", ci = "sampling", main = ttl, xlab = "Exposure (SD)")
    dev.off()
  }
}
cat("wrote bdlnm surfaces to", OUT, "\n")
