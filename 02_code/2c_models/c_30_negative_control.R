############ NEGATIVE CONTROL — DOES THE DESIGN ACTUALLY REMOVE SEASON? ############
# Earlier drafts of this diagnostic reported the spread of the climatological normal inside a
# stratum. That was a constructed quantity with no name in the literature, so it would have had
# to be defended in review. This does the standard thing instead (Lipsitch, Tchetgen Tchetgen &
# Cohen 2010): fit a NEGATIVE CONTROL EXPOSURE and see whether the design hands it an effect.
#
# The control is the 25-year climatological normal temperature for that cell and calendar week.
# It's pure season by construction - a 1997-2021 average contains no 2022 weather at all - and
# it is causally impossible, since a historical average cannot cause a 2022 outbreak. So there's
# no plausible-mechanism argument to have, which is rare for a negative control.
#
# Day length was the obvious alternative and is rejected: photoperiod drives waterfowl migration,
# so it has a real causal path to the outcome.
#
# A design that controls season returns ~1.00. The departure from 1 IS the residual seasonal
# confounding, as a rate ratio with an interval.
#
# Fitted with clogit rather than INLA on purpose: one linear term, no lag structure needed (the
# normal is smooth by construction), and conditional logistic is exact here. We already showed
# the conditional Poisson reproduces it to within 0.5% (c_28).

rm(list = ls())
project.folder = paste0(print(here::here()), '/')
source(paste0(project.folder, 'create_folder_structure.R'))
source(paste0(functions.folder, 'script_initiate.R'))
suppressMessages(library(survival))

# 1. the control itself: 25-year weekly normal temperature per cell
clim <- rbindlist(lapply(
  list.files(env_data, pattern = "flyway8_climatology_1997_2021_w", full.names = TRUE),
  fread, select = c("zone_id", "week_idx", "temperature_mean"),
  colClasses = list(character = "zone_id")))
setkey(clim, zone_id, week_idx)
cat(sprintf("climatology: %s cell-weeks, %d cells\n",
            format(nrow(clim), big.mark = ","), uniqueN(clim$zone_id)))

DESIGNS <- list(
  `Unidirectional, -7 to -28 d`                  = "case_crossover_df_multistate_daily.RDS",
  `Symmetric bidirectional, +/-14 to 28 d`       = "case_crossover_df_sym_g14.RDS",
  `Time-stratified, two-month strata`            = "case_crossover_df_timestrat_bimonth_nowo.RDS",
  `Time-stratified, one-month + 7-day post exclusion`   = "case_crossover_df_timestrat_month_post7.RDS")

one <- function(f, lab) {
  d <- setDT(readRDS(paste0(objects_folder, f)))
  if (anyDuplicated(names(d))) d <- d[, which(!duplicated(names(d))), with = FALSE]
  d[, date := as.Date(date)]
  d[, week_idx := pmin(as.integer(strftime(date, "%j")) %/% 7L, 51L)]
  d <- merge(d, clim, by = c("zone_id", "week_idx"), all.x = TRUE)
  d <- d[!is.na(temperature_mean)]

  # z-score within the panel so the contrast means the same thing across designs
  d[, nc_z := as.numeric(scale(temperature_mean))]
  # strata need a case and at least one referent to contribute anything
  d[, `:=`(k = sum(outbreak_binary), n = .N), by = stratum_id]
  d <- d[k == 1L & n > 1L]

  span <- d[, .(s = diff(range(as.integer(date)))), by = stratum_id][, mean(s)]
  fit <- clogit(outbreak_binary ~ nc_z + strata(stratum_id), data = d)
  b <- coef(fit)[["nc_z"]]; se <- sqrt(diag(vcov(fit)))[["nc_z"]]
  # reported per +0.5 SD, matching every other contrast in the paper
  data.table(Design = lab,
             Cases = sum(d$outbreak_binary),
             `Referents per case` = sprintf("%.1f", sum(d$outbreak_binary == 0) / sum(d$outbreak_binary)),
             `Stratum span (days)` = sprintf("%.1f", span),
             `Negative control RR` = sprintf("%.2f (%.2f, %.2f)",
                                             exp(0.5 * b), exp(0.5 * (b - 1.96 * se)),
                                             exp(0.5 * (b + 1.96 * se))),
             p = sprintf("%.3g", 2 * pnorm(-abs(b / se))))
}

out <- rbindlist(lapply(names(DESIGNS), function(nm) {
  r <- one(DESIGNS[[nm]], nm)
  cat(sprintf("  %-44s n=%3d  span %5s d  NC %s\n", nm, r$Cases,
              r$`Stratum span (days)`, r$`Negative control RR`))
  r
}))

saveRDS(out, paste0(objects_folder, "si_negative_control.RDS"))
cat("\n=== negative control: 25-year climatological normal temperature, per +0.5 SD ===\n")
print(out)
cat("\nA design that removes season should return ~1.00.\n")
