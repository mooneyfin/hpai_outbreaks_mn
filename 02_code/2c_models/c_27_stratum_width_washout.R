############ SUPPLEMENT — HOW WIDE SHOULD A STRATUM BE, AND HOW LONG THE WASHOUT? ############
# Two free knobs in the referent design, swept together because they trade off against each
# other. Neither is a "human" convention: dropping day-of-week matching is the only part of the
# usual time-stratified recipe that exists for human activity cycles, and we already dropped it.
# What's left is a fixed calendar partition, which is what keeps the conditional likelihood
# unbiased (Janes, Sheppard & Lumley 2005), and that argument doesn't care what the host is.
#
# Stratum width: wider blocks buy referents, but let more season back inside the stratum.
# Washout:       post-only, since only the days AFTER a detection leave the risk set.
#
# The post-exclusion length has an empirical anchor in these data. Within-cell repeat detections are
# bimodal: 9 of 14 land within 14 days (median 10.5), then nothing at all until 187 days. And
# the presumptive-to-confirmed reporting lag is median 2 days, 95th percentile 5. So anything
# from ~5 to 14 days sits in a defensible range and 14 clears the whole local-spread cluster.
#
# Reported per configuration: referents per case, how much 25-year seasonal signal survives
# inside a stratum, and the three terms the paper actually leans on.

rm(list = ls())
project.folder = paste0(print(here::here()), '/')
source(paste0(project.folder, 'create_folder_structure.R'))
source(paste0(functions.folder, 'script_initiate.R'))
suppressMessages({library(INLA); library(Matrix); library(flextable); library(officer)})
INLA::inla.setOption(num.threads = 1)
source(paste0(functions.folder, 'inla_dlnm_helpers.R'))

PREC <- 4; NS2 <- list(fun = "ns", df = 2); LIN <- list(fun = "lin")
SH <- c("temperature_shock", "precipitation_shock", "soil_moisture_shock", "wind_speed_shock")
AV <- c(setNames(rep(list(NS2), length(SH)), SH), list(anseriformes = LIN))
cell <- function(e) sprintf("%.2f (%.2f, %.2f)", e[["OR"]], e[["low"]], e[["high"]])

# 1. borrow the panel builder from a_06 rather than copy it. everything above its `for (dw in
#    ...)` export loop is setup plus build_timestrat() itself. two lines have to go first:
#    rm(list = ls()) and the create_folder_structure source, which ALSO starts with rm(list =
#    ls()) - and source() ignores the env you call it from, so it would wipe this script's
#    globals no matter where we eval. the folder vars are already set above anyway.
a06 <- readLines(paste0(data.prep.folder, "a_06_build_symmetric_panels.R"))
a06 <- a06[1:(grep("^for \\(dw in", a06)[1] - 1L)]
a06 <- a06[!grepl("^\\s*(rm\\(list|project\\.folder\\s*=|source\\()", a06)]
a06env <- new.env(parent = globalenv())
eval(parse(text = paste(a06, collapse = "\n")), envir = a06env)
build_timestrat <- get("build_timestrat", a06env)   # closure keeps a06env, so its helpers come too
daily <- get("daily", a06env)
stopifnot(is.function(build_timestrat), nrow(daily) > 0)
cat(sprintf("builder loaded | daily rows %s\n", format(nrow(daily), big.mark = ",")))

# 2. 25-year weekly normals, for the residual-season measure
clim <- rbindlist(lapply(
  list.files(env_data, pattern = "flyway8_climatology_1997_2021_w", full.names = TRUE),
  fread, select = c("zone_id", "week_idx", "temperature_mean")))
tot_t <- sd(clim$temperature_mean, na.rm = TRUE)

# 3. 90-day shock columns, attached onto whatever panel we build
ex <- setDT(readRDS(paste0(objects_folder, "daily_shock_predictors.RDS")))
ex[, date := as.Date(date)]
SHOCK <- setdiff(names(ex), c("zone_id", "date"))
attach_shocks <- function(cc) {
  if (anyDuplicated(names(cc))) cc <- cc[, which(!duplicated(names(cc))), with = FALSE]
  cc[, date := as.Date(date)]
  for (L in 0:28) {
    s <- ex[, c(.(zone_id = zone_id, date = date + L), .SD), .SDcols = SHOCK]
    setnames(s, SHOCK, paste0(SHOCK, "_Lag", L))
    cc <- merge(cc, s, by = c("zone_id", "date"), all.x = TRUE, sort = FALSE)
  }
  for (v in SHOCK) for (L in 0:28) {
    cn <- paste0(v, "_Lag", L); set(cc, j = cn, value = as.numeric(scale(cc[[cn]])))
  }
  need <- paste0(rep(SHOCK, each = 29), "_Lag", 0:28)
  cc <- cc[complete.cases(cc[, ..need])]
  cc[, `:=`(nc = sum(outbreak_binary), nr = sum(outbreak_binary == 0)), by = stratum_id]
  cc[nc == 1L & nr > 0L]
}

# 4. the sweep. washout 0 means every day in the block is a referent.
GRID <- CJ(bm = c(1L, 2L), wo = c(0L, 3L, 5L, 7L, 10L, 14L))[order(bm, wo)]

out <- rbindlist(lapply(seq_len(nrow(GRID)), function(i) {
  bm <- GRID$bm[i]; wo <- GRID$wo[i]
  raw <- build_timestrat(daily, block_months = bm, washout = wo, washout_type = "post")
  P <- attach_shocks(raw)

  # residual season: spread of the pure-seasonal 25-yr normal WITHIN a stratum
  z <- P[, .(zone_id, date, stratum_id)]
  z[, week_idx := pmin(as.integer(strftime(date, "%j")) %/% 7L, 51L)]
  m <- merge(z, clim, by = c("zone_id", "week_idx"), all.x = TRUE)
  s <- m[, .(span = diff(range(as.integer(date))),
             t_sd = sd(temperature_mean, na.rm = TRUE)), by = stratum_id]

  rbindlist(lapply(c("Minnesota", "Pooled"), function(pn) {
    d <- if (pn == "Minnesota") P[state == "Minnesota"] else copy(P)
    o <- assemble_cc_dlnm(d, met_vars = c(SH, "anseriformes"), weekly_vars = character(0),
                          max_lag = 28, argvar_by_var = AV, arglag = list(fun = "ns", df = 2))
    f <- suppressWarnings(fit_cc_inla_dlnm(o, fixed_prec = PREC))
    r <- data.table(
      `Stratum` = if (bm == 1L) "1 month" else "2 months",
      `Post exclusion (days)` = wo, Panel = pn,
      Cases = sum(d$outbreak_binary),
      `Referents per case` = sprintf("%.1f", sum(d$outbreak_binary == 0) / sum(d$outbreak_binary)),
      `Residual season (degC)` = sprintf("%.2f", mean(s$t_sd, na.rm = TRUE)),
      `% annual` = sprintf("%.0f%%", 100 * mean(s$t_sd, na.rm = TRUE) / tot_t),
      WAIC = sprintf("%.1f", f$waic$waic),
      `Temperature 0-7 days` = cell(met_lag_effect(f, o, "temperature_shock", 0:7)),
      `Soil moisture 0-7 days` = cell(met_lag_effect(f, o, "soil_moisture_shock", 0:7)),
      `Anseriformes 22-28 days` = cell(met_lag_effect(f, o, "anseriformes", 22:28)))
    cat(sprintf("  %-9s wo %2d | %-9s n=%3d  %5s ref/case  season %5s  WAIC %s\n",
                r$Stratum, wo, pn, r$Cases, r$`Referents per case`,
                r$`Residual season (degC)`, r$WAIC))
    r
  }))
}))

saveRDS(out, paste0(objects_folder, "tableS3_stratum_washout.RDS"))

is_sig <- function(x) {
  m <- regmatches(x, regexec("\\(([0-9.]+), ([0-9.]+)\\)", x))
  vapply(m, function(z) length(z) == 3 && (as.numeric(z[2]) > 1 | as.numeric(z[3]) < 1), TRUE)
}
ft <- flextable(out) |> fontsize(size = 9, part = "all") |> autofit() |>
  align(align = "center", part = "all") |> theme_booktabs()
for (j in grep("days$", names(out), value = TRUE)) {
  hit <- which(is_sig(out[[j]])); if (length(hit)) ft <- bold(ft, i = hit, j = j)
}
read_docx() |> body_add_flextable(ft) |>
  print(target = file.path(tables_main_folder, "TableS3_stratum_washout.docx"))

cat("\n=== stratum width x washout ===\n"); print(out, nrows = 100)
cat("\nwrote TableS3_stratum_washout.docx\n")
