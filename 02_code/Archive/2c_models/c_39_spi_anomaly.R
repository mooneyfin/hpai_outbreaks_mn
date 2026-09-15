############ GAMMA-BASED (SPI-STYLE) ANOMALY FOR THE SKEWED EXPOSURES ############
# The climatological anomaly is currently z = (x - mean) / sd, which assumes a roughly symmetric
# distribution. That holds for temperature (skew -0.24) and soil moisture (-0.59) and fails for
# precipitation (skew 4.09, max/median 505) and runoff (skew 13.7, 35% exact zeros). A z-score on
# runoff is close to meaningless.
#
# The standard fix in climate science is SPI (McKee et al. 1993): fit a gamma to the historical
# distribution, evaluate the observation's percentile under it, then map that percentile through
# qnorm so the result is on a familiar standard-normal scale.
#
# LIMITATION, and it is not small. The GEE export stores only _mean and _sd per cell-week, not
# the 25 individual years, so the gamma is fitted by METHOD OF MOMENTS rather than maximum
# likelihood, and the historical zero fraction is unknown. Proper SPI models zeros as a separate
# point mass; we cannot. So this is a genuine improvement on the z-score for precipitation and a
# partial one for runoff, not a textbook SPI. The GEE script should export percentiles directly
# if this is ever to be done properly.

rm(list = ls())
project.folder = paste0(print(here::here()), '/')
source(paste0(project.folder, 'create_folder_structure.R'))
source(paste0(functions.folder, 'script_initiate.R'))
suppressMessages({library(INLA); library(Matrix)})
INLA::inla.setOption(num.threads = 2)
source(paste0(functions.folder, 'inla_dlnm_helpers.R'))

SKEWED <- c("precipitation", "runoff")     # the two the z-score is wrong for
SYMM   <- c("temperature", "soil_moisture", "wind_speed")

clim <- rbindlist(lapply(
  list.files(env_data, pattern = "flyway8_climatology_1997_2021_w", full.names = TRUE),
  fread, colClasses = list(character = "zone_id")))
daily <- setDT(readRDS(paste0(objects_folder, multistate_ts_daily_rds)))
daily[, date := as.Date(date)][, week_idx := pmin(as.integer(strftime(date, "%j")) %/% 7L, 51L)]

# ERA5 depths are metres in both files, so no unit conversion is needed before comparing
d <- merge(daily[, c("zone_id", "date", "week_idx", SKEWED, SYMM), with = FALSE],
           clim[, c("zone_id", "week_idx",
                    paste0(c(SKEWED, SYMM), "_mean"), paste0(c(SKEWED, SYMM), "_sd")),
                with = FALSE],
           by = c("zone_id", "week_idx"), all.x = TRUE)

for (v in SKEWED) {
  m <- d[[paste0(v, "_mean")]]; s <- d[[paste0(v, "_sd")]]; x <- d[[v]]
  # gamma by method of moments: shape = (mean/sd)^2, rate = mean/sd^2
  shape <- (m / s)^2; rate <- m / s^2
  ok <- is.finite(shape) & is.finite(rate) & shape > 0 & rate > 0 & !is.na(x)
  pct <- rep(NA_real_, length(x))
  pct[ok] <- pgamma(pmax(x[ok], 0), shape = shape[ok], rate = rate[ok])
  pct <- pmin(pmax(pct, 1e-4), 1 - 1e-4)          # keep qnorm finite at the tails
  set(d, j = paste0(v, "_spi"), value = qnorm(pct))
}
for (v in SYMM) {                                   # unchanged: the z-score is appropriate here
  m <- d[[paste0(v, "_mean")]]; s <- d[[paste0(v, "_sd")]]
  set(d, j = paste0(v, "_spi"), value = fifelse(s > 0, (d[[v]] - m) / s, NA_real_))
}

cat("distribution of the new anomalies (should be ~N(0,1) if the fit is reasonable):\n")
for (v in c(SKEWED, SYMM)) {
  z <- d[[paste0(v, "_spi")]]
  cat(sprintf("  %-14s mean %+.2f  sd %.2f  skew %+.2f\n", v, mean(z, na.rm = TRUE),
              sd(z, na.rm = TRUE), mean(((z - mean(z, na.rm = TRUE)) /
                                         sd(z, na.rm = TRUE))^3, na.rm = TRUE)))
}
saveRDS(d[, c("zone_id", "date", paste0(c(SKEWED, SYMM), "_spi")), with = FALSE],
        paste0(objects_folder, "daily_spi_anomaly.RDS"))
cat("\nwrote daily_spi_anomaly.RDS\n")
