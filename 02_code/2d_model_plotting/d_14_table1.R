############ TABLE 1 — SPILLOVER CELL CHARACTERISTICS ############
# Three blocks: the premises in cells that had a spillover, land cover, and the exposures by
# meteorological season. The exposure block now reports BOTH parameterisations - the absolute
# condition and the anomaly - because the models use anomalies and a reader needs to see what
# an anomaly of a given size actually is in natural units.
#
# ERA5 ships temperature in KELVIN and depths in METRES; both are converted here or the table
# reports a mean temperature of 279 and a rainfall of 0.0002.

rm(list = ls())
project.folder = paste0(print(here::here()), '/')
source(paste0(project.folder, 'create_folder_structure.R'))
source(paste0(functions.folder, 'script_initiate.R'))
suppressMessages({library(flextable); library(officer)})

panel <- setDT(readRDS(paste0(objects_folder, multistate_ts_daily_rds)))[state == "Minnesota"]
scr   <- readRDS(paste0(objects_folder, multistate_independence_screening_rds))
count_flagged <- function(s, ids) {
  if (is.na(s)) return(0L)
  as.integer(sum(as.numeric(trimws(unlist(strsplit(s, ",")))) %in% ids))
}
panel[, nf := vapply(outbreak_id, count_flagged, integer(1), ids = scr$cross_farm_ids)]
panel[, outbreak_count := pmax(outbreak_count - nf, 0L)][, nf := NULL]
panel[, outbreak_binary := as.integer(outbreak_count > 0)]

case_cells <- unique(panel[outbreak_binary == 1]$zone_id)
cells <- unique(panel[zone_id %in% case_cells,
                      .(zone_id, feedlot_count, au_count, birds_au, cafo_flag)])
n_events <- sum(panel$outbreak_binary)

# 1. Premises
num_row <- function(x, lab, dp = 0) data.table(
  Characteristic = lab,
  Value = sprintf(paste0("%.", dp, "f (%.", dp, "f-%.", dp, "f)"),
                  median(x, na.rm = TRUE), quantile(x, .25, na.rm = TRUE),
                  quantile(x, .75, na.rm = TRUE)))
tab_premises <- rbindlist(list(
  data.table(Characteristic = "Spillover events", Value = as.character(n_events)),
  data.table(Characteristic = "Cells with at least one event", Value = as.character(length(case_cells))),
  num_row(cells$feedlot_count, "Registered feedlots per cell"),
  num_row(cells$au_count,      "Animal units per cell"),
  num_row(cells$birds_au,      "Birds per cell"),
  data.table(Characteristic = "Cells with a CAFO",
             Value = sprintf("%d (%.0f%%)", sum(cells$cafo_flag, na.rm = TRUE),
                             100 * mean(cells$cafo_flag, na.rm = TRUE)))))

# 2. Land cover. NLCD was only ever extracted on the legacy MN fishnet, so it comes across by
#    areal overlap between the two 10 km grids.
lcm <- setDT(readRDS(paste0(objects_folder, "mn8_land_cover_lookup.RDS")))
lc <- merge(cells[, .(zone_id)], lcm[, .(zone_id, land_cover_label)], by = "zone_id")[
  !is.na(land_cover_label), .N, by = land_cover_label][order(-N)]
tab_lc <- if (nrow(lc)) {
  lc[, Value := sprintf("%d (%.0f%%)", N, 100 * N / sum(N))]
  lc[, .(Characteristic = land_cover_label, Value)]
} else data.table(Characteristic = "Not available", Value = "-")

# 3. Events by meteorological season
ev <- panel[outbreak_binary == 1, .N, by = season][order(-N)]
ev[, Value := sprintf("%d (%.0f%%)", N, 100 * N / sum(N))]
tab_season <- ev[, .(Characteristic = as.character(season), Value)]

# 4. Exposures by season, absolute and anomaly, on the spillover cells
MET <- c(temperature = "Temperature (degC)", soil_moisture = "Soil moisture (m3/m3)",
         runoff = "Surface runoff (mm)", wind_speed = "Wind speed (m/s)")
LOGV <- c("runoff")
d <- setDT(readRDS(paste0(objects_folder, multistate_ts_daily_rds)))[
  state == "Minnesota" & zone_id %in% case_cells]
setorder(d, zone_id, date); d[, date := as.Date(date)]
d[, temperature := temperature - 273.15]
for (v in c("runoff", "precipitation")) set(d, j = v, value = d[[v]] * 1000)
for (v in names(MET)) {
  x <- if (v %in% LOGV) log1p(pmax(d[[v]], 0)) else d[[v]]
  set(d, j = "..x", value = x)
  d[, (paste0(v, "_anom")) := frollmean(..x, 7, align = "right") -
      shift(frollmean(..x, 90, align = "right"), 1), by = zone_id]
}
d[, ..x := NULL]

fmt <- function(x, dp) sprintf(paste0("%.", dp, "f (%.", dp, "f, %.", dp, "f)"),
                               mean(x, na.rm = TRUE), quantile(x, .05, na.rm = TRUE),
                               quantile(x, .95, na.rm = TRUE))
rows <- rbindlist(lapply(names(MET), function(v) {
  dp <- if (v == "soil_moisture") 3 else 1
  a <- d[, .(V = fmt(get(v), dp)), by = season]
  b <- d[, .(V = fmt(get(paste0(v, "_anom")), if (v %in% LOGV) 2 else dp)), by = season]
  rbind(cbind(data.table(Characteristic = paste0(MET[[v]], " - absolute")), a),
        cbind(data.table(Characteristic = paste0(MET[[v]], " - anomaly")),  b))
}))
tab_met <- dcast(rows, Characteristic ~ season, value.var = "V")
setcolorder(tab_met, c("Characteristic",
                       intersect(c("Winter", "Spring", "Summer", "Fall"), names(tab_met))))

# 5. Assemble
sec <- function(t, lab) rbindlist(list(data.table(Characteristic = lab, Value = ""), t), fill = TRUE)
tab1 <- rbindlist(list(sec(tab_premises, "Premises"), sec(tab_lc, "Land cover"),
                       sec(tab_season, "Events by season")), fill = TRUE)
ft <- flextable(tab1) |> bold(i = which(tab1$Value == ""), j = 1) |>
  fontsize(size = 11, part = "all") |> autofit() |> theme_booktabs() |>
  align(j = 2, align = "right", part = "all")
ft_met <- flextable(tab_met) |> fontsize(size = 11, part = "all") |> autofit() |>
  theme_booktabs() |> align(j = 2:ncol(tab_met), align = "center", part = "all")

saveRDS(list(cells = tab1, exposures = tab_met), paste0(objects_folder, "table1_parts.RDS"))
read_docx() |> body_add_flextable(ft) |>
  print(target = file.path(tables_main_folder, "Table1_spillover_cells.docx"))
read_docx() |> body_add_flextable(ft_met) |>
  print(target = file.path(tables_main_folder, "Table2_exposures_by_season.docx"))

print(tab1); cat("\n"); print(tab_met)
cat("\nwrote Table1_spillover_cells.docx and Table2_exposures_by_season.docx\n")
