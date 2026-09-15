############ DOES THE WATERFOWL FINDING NEED THE LOG TRANSFORM? ############
# Bird abundance and precipitation are right-skewed, so the panel applies log1p before
# z-scoring (a_06, `skewed_vars`). On the SUPERSEDED specification that choice moved the
# headline: waterfowl read 1.66 (1.23, 2.26) logged against 1.02 (0.74, 1.40) raw. That check
# was never repeated on the current spec, and "the result requires log1p" is the first thing a
# reviewer will ask.
#
# The panel only stores the transformed columns, so the raw-scale versions get rebuilt here
# from the daily series and swapped in. Everything else - design, washout, prior, crossbasis,
# lag window - is held fixed, so the only thing that changes between rows is the transform.

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
MAXLAG <- 28
cell <- function(e) sprintf("%.2f (%.2f, %.2f)", e[["OR"]], e[["low"]], e[["high"]])

D <- setDT(readRDS(paste0(objects_folder, "case_crossover_df_shock_post7.RDS")))
if (anyDuplicated(names(D))) D <- D[, which(!duplicated(names(D))), with = FALSE]
D[, date := as.Date(date)]

# 1. Rebuild the bird lag columns on the RAW scale. Same lag construction and same z-scoring as
#    the panel uses, just without the log1p, so the two differ in exactly one step.
daily <- setDT(readRDS(paste0(objects_folder, multistate_ts_daily_rds)))
setorder(daily, zone_id, date)
raw <- daily[, .(zone_id, date = as.Date(date), anseriformes)]
nm  <- paste0("rawbird_Lag", 0:MAXLAG)
raw[, (nm) := shift(anseriformes, n = 0:MAXLAG, type = "lag"), by = zone_id]
D <- merge(D, raw[, c("zone_id", "date", nm), with = FALSE],
           by = c("zone_id", "date"), all.x = TRUE, sort = FALSE)
for (cn in nm) set(D, j = cn, value = as.numeric(scale(D[[cn]])))

# keep only rows where BOTH versions exist, so logged and raw are fitted on identical data
D <- D[complete.cases(D[, ..nm])]
D[, `:=`(nc = sum(outbreak_binary), nr = sum(outbreak_binary == 0)), by = stratum_id]
D <- D[nc == 1L & nr > 0L]
cat(sprintf("common rows: %d cases, %d strata\n", sum(D$outbreak_binary), uniqueN(D$stratum_id)))

# rawbird_Lag* has to be visible to the crossbasis builder under a variable name, and
# assemble_cc_dlnm looks for "<var>_Lag<k>", so the prefix IS the variable name
fit_one <- function(dat, birdvar) {
  mv <- c(SH, birdvar)
  av <- c(setNames(rep(list(NS2), length(SH)), SH), setNames(list(LIN), birdvar))
  o  <- assemble_cc_dlnm(dat, met_vars = mv, weekly_vars = character(0), max_lag = MAXLAG,
                         argvar_by_var = av, arglag = list(fun = "ns", df = 2))
  f  <- suppressWarnings(fit_cc_inla_dlnm(o, fixed_prec = PREC))
  list(waic = f$waic$waic,
       b0 = cell(met_lag_effect(f, o, birdvar, 0:7)),
       b3 = cell(met_lag_effect(f, o, birdvar, 22:28)),
       sm = cell(met_lag_effect(f, o, "soil_moisture_shock", 0:7)))
}

PANELS <- c("Minnesota", "Other northern Mississippi Flyway states", "Pooled")
sub <- function(pn) {
  if (pn == "Minnesota") D[state == "Minnesota"]
  else if (pn == "Pooled") copy(D)
  else D[state != "Minnesota"]
}
VARIANTS <- c(`log1p, then z-scored (primary)` = BIRD, `raw counts, z-scored only` = "rawbird")

out <- rbindlist(lapply(PANELS, function(pn) {
  d <- sub(pn)
  rbindlist(lapply(names(VARIANTS), function(vn) {
    r <- fit_one(d, VARIANTS[[vn]])
    cat(sprintf("  %-42s %-30s WAIC %.1f  bird 22-28 %s\n", pn, vn, r$waic, r$b3))
    data.table(Panel = pn, `Bird abundance scale` = vn, Cases = sum(d$outbreak_binary),
               WAIC = sprintf("%.1f", r$waic),
               `Anseriformes 0-7 days` = r$b0, `Anseriformes 22-28 days` = r$b3,
               `Soil moisture 0-7 days` = r$sm)
  }))
}))

saveRDS(out, paste0(objects_folder, "tableS9_log_transform.RDS"))
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
  print(target = file.path(tables_main_folder, "TableS9.docx"))

print(out)
cat("\nwrote TableS9.docx\n")
