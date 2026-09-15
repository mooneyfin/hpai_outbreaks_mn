# a_06_build_symmetric_panels.R
# Symmetric bidirectional case-crossover panels, lags out to 56 days.
#
# Why this exists: the symmetric panels c_11 fits were built in a scratch script that never
# made it into the repo, and they only carried lags 0-28. Two things needed fixing:
#   - the original manuscript ran the bird term out to 8 lag weeks and found the waterfowl
#     peak around week 5-6. At a 28-day cap our curve is still climbing at the boundary, so
#     the window was probably truncating the effect.
#   - the +/-14 d referent gap was never compared against +/-7 on the current panel.
#
# The panels STORE 56 lag columns; every model still fits max_lag = 28. The deeper columns
# exist only so the 42/56 lag-window sensitivity can be fitted without a rebuild. Because each
# lag column is z-scaled on its own, columns 0-28 are identical whether the file holds 28 or
# 56 of them, which is why these reproduce the previous panels exactly.

# 0a. Root, folders, packages
rm(list = ls())
project.folder = paste0(print(here::here()), '/')
source(paste0(project.folder, 'create_folder_structure.R'))
source(paste0(functions.folder, 'script_initiate.R'))

LAG_STORE <- 56   # columns stored; models fit max_lag = 28
FAR      <- 28                     # outermost referent, days either side
GAPS     <- c(7, 14)               # innermost referent (the washout)

# 1a. Inputs
daily     <- setDT(readRDS(paste0(objects_folder, "multistate_fishnet_timeseries_daily.RDS")))
screening <- readRDS(paste0(objects_folder, "multistate_independence_screening.RDS"))
setkey(daily, zone_id, date)

case_window_start <- as.Date("2022-01-01")
case_window_end   <- as.Date("2022-12-31")
cross_farm_ids    <- screening$cross_farm_ids

count_flagged_ids <- function(id_str, flagged_ids) {
  if (is.na(id_str)) return(0L)
  as.integer(sum(as.numeric(trimws(unlist(strsplit(id_str, ",")))) %in% flagged_ids))
}

lag_vars <- intersect(c(
  "temperature", "dewpoint", "relative_humidity", "wind_speed",
  "skin_temperature", "soil_temperature", "soil_moisture", "soil_moisture_2",
  "snow_depth", "snow_cover", "snow_density",
  "precipitation", "snowfall", "runoff", "evapotranspiration",
  "leaf_area_index", "shortwave_radiation", "longwave_radiation", "surface_pressure",
  "total_bird_abundance", "anseriformes", "predators"), names(daily))
skewed_vars <- intersect(c("total_bird_abundance", "anseriformes", "predators",
                           "runoff", "precipitation"), lag_vars)

# ERA5 ships runoff and precipitation in METRES, order 1e-4, where log1p is essentially the
# identity and corrects none of the skew. The original panels converted to mm first, and
# skipping that step silently changes the exposure scale: max z goes 11.4 -> 32.0 for runoff.
DEPTH_M_TO_MM <- intersect(c("runoff", "precipitation", "snowfall"), lag_vars)

transform_lag_cols <- function(dt, max_lag, use_log = TRUE) {
  for (v in DEPTH_M_TO_MM) for (k in 0:max_lag) {
    col <- paste0(v, "_Lag", k); if (col %in% names(dt)) dt[[col]] <- dt[[col]] * 1000
  }
  if (use_log) for (v in skewed_vars) for (k in 0:max_lag) {
    col <- paste0(v, "_Lag", k); if (col %in% names(dt)) dt[[col]] <- log1p(dt[[col]])
  }
  for (v in lag_vars) for (k in 0:max_lag) {
    col <- paste0(v, "_Lag", k); if (col %in% names(dt)) dt[[col]] <- as.numeric(scale(dt[[col]]))
  }
  dt
}

# 2a. Same shape as build_uni_cc in a_03, but referents sit on BOTH sides of the case at
#     +/-gap .. +/-far, and the lag matrix runs deeper.
build_sym_cc <- function(panel, gap, far = FAR, max_lag = LAG_STORE,
                         flag_ids = cross_farm_ids, use_log = TRUE) {
  p <- copy(panel)
  p[, tkey := as.Date(date)]
  p[, is_outbreak_raw := as.integer(outbreak_binary == 1)]

  # case definition: drop cross-farm events, then redefine cases on what's left
  fl <- p[outbreak_binary == 1 & tkey >= case_window_start & tkey <= case_window_end
        ][, nf := vapply(outbreak_id, count_flagged_ids, integer(1), flagged_ids = flag_ids)
        ][nf > 0, .(zone_id, tkey, nf)]
  p <- merge(p, fl, by = c("zone_id", "tkey"), all.x = TRUE)
  p[, nf := fcoalesce(nf, 0L)]
  p[, outbreak_binary := as.integer(outbreak_count - nf > 0)]
  p[, nf := NULL]

  cases <- p[outbreak_binary == 1 & tkey >= case_window_start & tkey <= case_window_end,
             .(zone_id, ct = tkey)]
  cases[, stratum_id := paste0("z", zone_id, "_", format(ct, "%Y%m%d"))]

  # lag matrix on the case zones' full series
  src <- p[zone_id %in% unique(cases$zone_id)]
  setorder(src, zone_id, tkey)
  for (v in lag_vars) {
    nm <- paste0(v, "_Lag", 0:max_lag)
    src[, (nm) := shift(.SD[[1]], n = 0:max_lag, type = "lag"), by = zone_id, .SDcols = v]
  }

  offs <- c(0L, -(gap:far), (gap:far))          # case day plus both referent arms
  cases[, cid := .I]
  sets <- cases[, .(offset = offs, tkey = ct + offs,
                    is_case = c(1L, rep(0L, length(offs) - 1L))),
                by = .(cid, zone_id, stratum_id)]

  matched <- merge(sets, src, by = c("zone_id", "tkey"), all.x = FALSE)
  # a referent day that was itself an outbreak isn't a valid control
  matched <- matched[is_case == 1L | is_outbreak_raw == 0L]
  keep <- matched[, .(nc = sum(is_case == 1L), nk = sum(is_case == 0L)),
                  by = stratum_id][nc >= 1 & nk >= 1, stratum_id]
  matched <- matched[stratum_id %in% keep]

  matched[, `:=`(outbreak_binary = is_case,
                 case_control    = fifelse(is_case == 1L, "case", "control"))]
  setnames(matched, "tkey", "date")
  transform_lag_cols(matched, max_lag, use_log)
}

# 3a. Build and save
for (g in GAPS) {
  d <- build_sym_cc(daily, gap = g)
  f <- sprintf("case_crossover_df_sym_g%d.RDS", g)
  saveRDS(d, paste0(objects_folder, f))
  mn <- d[state == "Minnesota"]
  cat(sprintf("gap +/-%2d: %6d rows | MN cases %3d, referents %5d (%.1f per case), cells %3d -> %s\n",
              g, nrow(d), sum(mn$outbreak_binary == 1), sum(mn$outbreak_binary == 0),
              sum(mn$outbreak_binary == 0) / sum(mn$outbreak_binary == 1),
              uniqueN(mn$zone_id), f))
}

# 3b. Sensitivity variants, all at the primary gap so only one thing changes at a time.
#     - no-log: is log1p doing anything scaling alone doesn't (note m vs mm is irrelevant
#       without the log, since z-scaling is invariant to a constant multiplier)
#     - CSLT-inclusive: keep the cross-farm events the primary drops
#     - 200 m / 1 km: add the spatio-temporal independence screens on top of CSLT
VARIANTS <- list(
  nolog = list(flag = cross_farm_ids,                              use_log = FALSE),
  cslt  = list(flag = integer(0),                                  use_log = TRUE),
  s200m = list(flag = union(cross_farm_ids, screening$ids_200m),   use_log = TRUE),
  s1km  = list(flag = union(cross_farm_ids, screening$ids_1km),    use_log = TRUE)
)
for (nm in names(VARIANTS)) {
  v <- VARIANTS[[nm]]
  d <- build_sym_cc(daily, gap = 14, flag_ids = v$flag, use_log = v$use_log)
  f <- sprintf("case_crossover_df_sym_g14_%s.RDS", nm)
  saveRDS(d, paste0(objects_folder, f))
  mn <- d[state == "Minnesota"]
  cat(sprintf("%-6s MN cases %3d, referents %5d, cells %3d -> %s\n",
              nm, sum(mn$outbreak_binary == 1), sum(mn$outbreak_binary == 0),
              uniqueN(mn$zone_id), f))
}

# 4a. Weekly-resolution panel for SI Table S4, the "did weekly -> daily change anything" check.
#     Same design one step coarser: referents at +/-2 to +/-4 weeks (the week equivalent of
#     +/-14 to +/-28 days), lags 0-4 weeks (0-28 days).
weekly <- setDT(readRDS(paste0(objects_folder, "multistate_fishnet_timeseries_weekly.RDS")))
setkey(weekly, zone_id, week_start)

build_sym_weekly <- function(panel, gap = 2, far = 4, max_lag = 8,
                             flag_ids = cross_farm_ids) {
  p <- copy(panel)
  p[, tkey := as.Date(week_start)]
  p[, is_outbreak_raw := as.integer(outbreak_binary == 1)]

  fl <- p[outbreak_binary == 1 & tkey >= case_window_start & tkey <= case_window_end
        ][, nf := vapply(outbreak_id, count_flagged_ids, integer(1), flagged_ids = flag_ids)
        ][nf > 0, .(zone_id, tkey, nf)]
  p <- merge(p, fl, by = c("zone_id", "tkey"), all.x = TRUE)
  p[, nf := fcoalesce(nf, 0L)]
  p[, outbreak_binary := as.integer(outbreak_count - nf > 0)]
  p[, nf := NULL]

  cases <- p[outbreak_binary == 1 & tkey >= case_window_start & tkey <= case_window_end,
             .(zone_id, ct = tkey)]
  cases[, stratum_id := paste0("z", zone_id, "_", format(ct, "%Y%m%d"))]

  src <- p[zone_id %in% unique(cases$zone_id)]
  setorder(src, zone_id, tkey)
  for (v in lag_vars) {
    nm <- paste0(v, "_Lag", 0:max_lag)
    src[, (nm) := shift(.SD[[1]], n = 0:max_lag, type = "lag"), by = zone_id, .SDcols = v]
  }

  offs <- c(0L, -(gap:far), (gap:far))
  cases[, cid := .I]
  sets <- cases[, .(tkey = ct + offs * 7L,                       # weeks -> days
                    is_case = c(1L, rep(0L, length(offs) - 1L))),
                by = .(cid, zone_id, stratum_id)]
  matched <- merge(sets, src, by = c("zone_id", "tkey"), all.x = FALSE)
  matched <- matched[is_case == 1L | is_outbreak_raw == 0L]
  keep <- matched[, .(nc = sum(is_case == 1L), nk = sum(is_case == 0L)),
                  by = stratum_id][nc >= 1 & nk >= 1, stratum_id]
  matched <- matched[stratum_id %in% keep]
  matched[, `:=`(outbreak_binary = is_case,
                 case_control    = fifelse(is_case == 1L, "case", "control"))]
  setnames(matched, "tkey", "week_start")
  transform_lag_cols(matched, max_lag)
}

wk <- build_sym_weekly(weekly)
saveRDS(wk, paste0(objects_folder, "case_crossover_df_sym_weekly.RDS"))
mnw <- wk[state == "Minnesota"]
cat(sprintf("weekly  MN cases %3d, referents %5d, cells %3d -> case_crossover_df_sym_weekly.RDS\n",
            sum(mnw$outbreak_binary == 1), sum(mnw$outbreak_binary == 0), uniqueN(mnw$zone_id)))


# 5a. Time-stratified panel, calendar-month strata, ALL days as referents.
#     The usual time-stratified recipe also matches day-of-week, which exists to control weekly
#     human activity cycles in air-pollution studies. That is not a mechanism here, and matching
#     on it costs most of the referents (3-4 per case instead of ~29) and with them the
#     within-stratum lag variation the bird term needs. Dropping it keeps the property that
#     matters - strata are a fixed calendar partition, so referent windows never overlap between
#     cases and the conditional likelihood stays unbiased (Janes, Sheppard & Lumley 2005).
build_timestrat <- function(panel, max_lag = LAG_STORE, flag_ids = cross_farm_ids,
                            match_dow = FALSE, block_months = 1L, washout = 0L,
                            washout_type = c("symmetric", "post")) {
  washout_type <- match.arg(washout_type)
  p <- copy(panel)
  p[, tkey := as.Date(date)]
  p[, is_outbreak_raw := as.integer(outbreak_binary == 1)]

  fl <- p[outbreak_binary == 1 & tkey >= case_window_start & tkey <= case_window_end
        ][, nf := vapply(outbreak_id, count_flagged_ids, integer(1), flagged_ids = flag_ids)
        ][nf > 0, .(zone_id, tkey, nf)]
  p <- merge(p, fl, by = c("zone_id", "tkey"), all.x = TRUE)
  p[, nf := fcoalesce(nf, 0L)]
  p[, outbreak_binary := as.integer(outbreak_count - nf > 0)]
  p[, nf := NULL]

  # the fixed calendar partition. block_months = 2 pairs calendar months (Jan-Feb, Mar-Apr,
  # ...), which keeps the partition fixed while giving a wider referent pool - useful once a
  # washout removes the days nearest the case.
  p[, blk := paste0(zone_id, "_", format(tkey, "%Y"), "_",
                    ceiling(as.integer(format(tkey, "%m")) / block_months),
                    if (match_dow) paste0("_", format(tkey, "%u")) else "")]

  cases <- p[outbreak_binary == 1 & tkey >= case_window_start & tkey <= case_window_end,
             .(zone_id, ct = tkey, blk)]
  cases[, stratum_id := paste0("z", zone_id, "_", format(ct, "%Y%m%d"))]

  src <- p[zone_id %in% unique(cases$zone_id)]
  setorder(src, zone_id, tkey)
  for (v in lag_vars) {
    nm <- paste0(v, "_Lag", 0:max_lag)
    src[, (nm) := shift(.SD[[1]], n = 0:max_lag, type = "lag"), by = zone_id, .SDcols = v]
  }

  # every day in the case's own block becomes a referent
  sets <- merge(cases[, .(stratum_id, zone_id, ct, blk)],
                src[, .(zone_id, tkey, blk)], by = c("zone_id", "blk"),
                allow.cartesian = TRUE)
  sets[, is_case := as.integer(tkey == ct)]
  matched <- merge(sets[, .(stratum_id, zone_id, tkey, is_case)], src,
                   by = c("zone_id", "tkey"), all.x = FALSE)
    # washout: referents within `washout` days of the case share most of their 0-28 d exposure
  # history with it, and after detection the flock is depopulated so those days are not valid
  # person-time either. Applied symmetrically to keep the referents balanced in time.
  if (washout > 0L) {
    matched <- merge(matched, cases[, .(stratum_id, ct)], by = "stratum_id", all.x = TRUE)
    dd <- as.integer(matched$tkey - matched$ct)
    matched <- matched[is_case == 1L |
                       if (washout_type == "symmetric") abs(dd) > washout else !(dd > 0 & dd <= washout)]
    matched[, ct := NULL]
  }
  matched <- matched[is_case == 1L | is_outbreak_raw == 0L]
  keep <- matched[, .(nc = sum(is_case == 1L), nk = sum(is_case == 0L)),
                  by = stratum_id][nc >= 1 & nk >= 1, stratum_id]
  matched <- matched[stratum_id %in% keep]
  matched[, `:=`(outbreak_binary = is_case,
                 case_control    = fifelse(is_case == 1L, "case", "control"))]
  setnames(matched, "tkey", "date")
  transform_lag_cols(matched, max_lag)
}

for (dw in c(FALSE, TRUE)) {
  d <- build_timestrat(daily, match_dow = dw)
  f <- sprintf("case_crossover_df_timestrat_%s.RDS", if (dw) "month_dow" else "month")
  saveRDS(d, paste0(objects_folder, f))
  mn <- d[state == "Minnesota"]
  cat(sprintf("time-stratified %-10s MN cases %3d, referents %5d (%.1f per case), cells %3d -> %s\n",
              if (dw) "month+DOW" else "month", sum(mn$outbreak_binary == 1),
              sum(mn$outbreak_binary == 0),
              sum(mn$outbreak_binary == 0) / sum(mn$outbreak_binary == 1),
              uniqueN(mn$zone_id), f))
}


# 5b. Two-month strata with a symmetric washout, built for the whole flyway panel so Minnesota
#     and the pooled analysis share one construction.
for (cfg in list(list(bm = 1L, wo = 14L, tag = "month_wo14"),
                 list(bm = 2L, wo = 14L, tag = "bimonth_wo14"),
                 list(bm = 2L, wo = 21L, tag = "bimonth_wo21"))) {
  d <- build_timestrat(daily, block_months = cfg$bm, washout = cfg$wo)
  f <- sprintf("case_crossover_df_timestrat_%s.RDS", cfg$tag)
  saveRDS(d, paste0(objects_folder, f))
  mn <- d[state == "Minnesota"]
  cat(sprintf("%-14s ALL: cases %3d refs %6d (%.1f/case) | MN: cases %3d refs %5d (%.1f/case)\n",
              cfg$tag, sum(d$outbreak_binary==1), sum(d$outbreak_binary==0),
              sum(d$outbreak_binary==0)/sum(d$outbreak_binary==1),
              sum(mn$outbreak_binary==1), sum(mn$outbreak_binary==0),
              sum(mn$outbreak_binary==0)/sum(mn$outbreak_binary==1)))
  print(d[outbreak_binary==1, .(cases=.N), by=state][order(-cases)])
}


# 5c. Two variants that separate the two competing explanations for the flattened lag profile:
#     is it exposure-history overlap between case and referent, or simply that a 28-day lag
#     window has no room inside a 30-day stratum?
for (cfg in list(list(bm = 2L, wo = 0L,  ty = "symmetric", tag = "bimonth_nowo"),
                 list(bm = 2L, wo = 14L, ty = "post",      tag = "bimonth_post14"))) {
  d <- build_timestrat(daily, block_months = cfg$bm, washout = cfg$wo, washout_type = cfg$ty)
  saveRDS(d, paste0(objects_folder, sprintf("case_crossover_df_timestrat_%s.RDS", cfg$tag)))
  mn <- d[state == "Minnesota"]
  cat(sprintf("%-16s ALL cases %3d | MN cases %3d refs %5d (%.1f/case)\n", cfg$tag,
              sum(d$outbreak_binary == 1), sum(mn$outbreak_binary == 1),
              sum(mn$outbreak_binary == 0),
              sum(mn$outbreak_binary == 0) / sum(mn$outbreak_binary == 1)))
}


# 5c-bis. Post-only washouts at one-month strata. The risk-set argument only justifies
#     dropping days AFTER the case, since the flock is depopulated and can no longer produce
#     the outcome; days before the case are valid person-time and there's no validity reason
#     to bin them. these panels let us check what the symmetric version actually buys.
for (cfg in list(list(bm = 1L, wo = 7L,   tag = "month_post7"),
                 list(bm = 1L, wo = 14L,  tag = "month_post14"),
                 list(bm = 1L, wo = 400L, tag = "month_postall"),
                 list(bm = 2L, wo = 400L, tag = "bimonth_postall"))) {
  d <- build_timestrat(daily, block_months = cfg$bm, washout = cfg$wo, washout_type = "post")
  saveRDS(d, paste0(objects_folder, sprintf("case_crossover_df_timestrat_%s.RDS", cfg$tag)))
  mn <- d[state == "Minnesota"]
  cat(sprintf("%-14s ALL cases %3d | MN cases %3d refs %5d (%.1f/case)\n", cfg$tag,
              sum(d$outbreak_binary == 1), sum(mn$outbreak_binary == 1),
              sum(mn$outbreak_binary == 0),
              sum(mn$outbreak_binary == 0) / sum(mn$outbreak_binary == 1)))
}


# 5d. Shorter washouts. 7 days clears depopulation while keeping far more referents than 14,
#     so it tests whether the de-attenuation needs a wide separation or just a few days.
for (cfg in list(list(bm = 1L, wo = 7L,  tag = "month_wo7"),
                 list(bm = 2L, wo = 7L,  tag = "bimonth_wo7"),
                 list(bm = 2L, wo = 10L, tag = "bimonth_wo10"))) {
  d <- build_timestrat(daily, block_months = cfg$bm, washout = cfg$wo)
  saveRDS(d, paste0(objects_folder, sprintf("case_crossover_df_timestrat_%s.RDS", cfg$tag)))
  mn <- d[state == "Minnesota"]
  cat(sprintf("%-14s ALL cases %3d | MN cases %3d refs %5d (%.1f/case)\n", cfg$tag,
              sum(d$outbreak_binary == 1), sum(mn$outbreak_binary == 1),
              sum(mn$outbreak_binary == 0),
              sum(mn$outbreak_binary == 0) / sum(mn$outbreak_binary == 1)))
}
