# d_08_consort_flowchart.R
# Sample-size cascade for the daily case-crossover study of HPAI spillover into MN poultry.
#
# Replaces d_04_consort_flowchart.R, which described the time-stratified design
# (zone x year x month x day-of-week) on the legacy MN fishnet. Both are superseded: the
# primary is now symmetric bidirectional referents at +/-14 to +/-28 d on the 8-state grid,
# and the old cascade objects say 88 events / 70 zones against the current 108 -> 85 / 73.
#
# Everything below is counted from the panels themselves, so it can't drift out of sync with
# the models the way the hardcoded cascade objects did.

# 0a. Root, folders, packages
rm(list = ls())
project.folder = paste0(print(here::here()), '/')
source(paste0(project.folder, 'create_folder_structure.R'))
source(paste0(functions.folder, 'script_initiate.R'))
source(paste0(functions.folder, 'inla_dlnm_helpers.R'))

# 1a. Raw detections on the analysis grid
panel <- setDT(readRDS(paste0(objects_folder, multistate_ts_daily_rds)))[state == "Minnesota"]
scr   <- readRDS(paste0(objects_folder, multistate_independence_screening_rds))

n_raw_events    <- sum(panel$outbreak_count, na.rm = TRUE)
n_raw_celldays  <- sum(panel$outbreak_count > 0, na.rm = TRUE)
n_raw_cells     <- uniqueN(panel[outbreak_count > 0]$zone_id)
n_multi_cellday <- sum(panel$outbreak_count > 1, na.rm = TRUE)
date_from       <- format(min(panel$date), "%d %b %Y")
date_to         <- format(max(panel$date), "%d %b %Y")

# 1b. Independence screen: drop events flagged as farm-to-farm spread rather than
#     independent introductions from the environment
count_flagged <- function(s, ids) {
  if (is.na(s)) return(0L)
  as.integer(sum(as.numeric(trimws(unlist(strsplit(s, ",")))) %in% ids))
}
panel[, nf := vapply(outbreak_id, count_flagged, integer(1), ids = scr$cross_farm_ids)]
n_flagged <- sum(panel$nf)
panel[, oc := pmax(outbreak_count - nf, 0L)]

n_events   <- sum(panel$oc)
n_celldays <- sum(panel$oc > 0)
n_cells    <- uniqueN(panel[oc > 0]$zone_id)

# 1c. Final analytic panel
cc <- setDT(readRDS(paste0(objects_folder, "case_crossover_df_sym_g14.RDS")))[
  state == "Minnesota"]
n_cases   <- sum(cc$outbreak_binary == 1)
n_refs    <- sum(cc$outbreak_binary == 0)
n_strata  <- uniqueN(cc$stratum_id)
n_cc_cells<- uniqueN(cc$zone_id)
per_case  <- n_refs / n_cases

# referent window is -28..-14 and +14..+28 d, so 15 days either side
max_refs      <- 2 * 15
n_refs_lost   <- n_cases * max_refs - n_refs
fmt <- function(x) formatC(x, big.mark = ",", format = "d")

# 2a. Cascade boxes, top to bottom
boxes <- data.table(
  # evenly spaced, boxes tall enough that the gaps between them don't dominate the figure
  y = seq(0.93, 0.07, length.out = 6),
  h = 0.115,
  label = c(
    sprintf("USDA APHIS confirmed HPAI poultry detections\nMinnesota, %s to %s\nn = %d events",
            date_from, date_to, n_raw_events),
    sprintf("Aggregated to 10 km × 1 day grid cells\n%d events → %d case cell-days in %d cells",
            n_raw_events, n_raw_celldays, n_raw_cells),
    sprintf("Independence screen\nCross-farm transmission events excluded\n− %d events",
            n_flagged),
    sprintf("Eligible cases\n%d events → %d case cell-days in %d cells",
            n_events, n_celldays, n_cells),
    sprintf("Referent selection within the same cell\nDays −28 to −14 and +14 to +28 relative to the case\nup to %d referents per case",
            max_refs),
    sprintf("Analytic sample\n%d case cell-days  |  %s referent cell-days\n%d strata across %d cells (%.1f referents per case)",
            n_cases, fmt(n_refs), n_strata, n_cc_cells, per_case))
)

# 2b. Arrows between consecutive boxes
arrows <- data.table(
  y    = boxes$y[-nrow(boxes)] - boxes$h[-nrow(boxes)] / 2,
  yend = boxes$y[-1]           + boxes$h[-1] / 2
)

# 2c. Right-hand column explaining what each step does. Each note sits level with the box it
#     explains; the last box needs none.
sides <- data.table(
  y = boxes$y[1:5],
  txt = c(
    "Detections with a confirmed outbreak date,\nmatched to the 10 km analysis grid",
    sprintf("%d cell-days carried more than one\ndetection and collapse to a single case", n_multi_cellday),
    "Flagged by sequence and epidemiological\ninvestigation as farm-to-farm spread,\nso not independent introductions",
    "The primary case definition.\nCSLT-inclusive counts are reported\nas a sensitivity analysis",
    sprintf("Symmetry cancels a linear time trend;\nthe ±14 d gap clears the depopulation\nwindow. %d referents fall outside the\nstudy period and are dropped", n_refs_lost))
)

# 3a. Draw it. No title or caption by design; those live in the manuscript text.
fp <- ggplot() +
  geom_rect(data = boxes,
            aes(xmin = 0.03, xmax = 0.60, ymin = y - h / 2, ymax = y + h / 2),
            fill = "white", colour = "grey25", linewidth = 0.4) +
  geom_text(data = boxes, aes(x = 0.315, y = y, label = label),
            family = "Avenir", size = 3.1, lineheight = 1.12) +
  geom_segment(data = arrows, aes(x = 0.315, xend = 0.315, y = y, yend = yend),
               arrow = arrow(length = unit(0.18, "cm"), type = "closed"),
               colour = "grey25", linewidth = 0.4) +
  geom_segment(aes(x = 0.65, xend = 0.65, y = 0.02, yend = 0.98),
               colour = "grey80", linewidth = 0.3, linetype = "dotted") +
  geom_text(data = sides, aes(x = 0.68, y = y, label = txt),
            family = "Avenir", size = 2.6, hjust = 0, colour = "grey30", lineheight = 1.12) +
  scale_x_continuous(limits = c(0, 1), expand = c(0, 0)) +
  scale_y_continuous(limits = c(0, 1), expand = c(0, 0)) +
  theme_spark_map()

# 4a. Export. PDF too, since journals prefer vector for line art — but through cairo_pdf, as
#     the default pdf() device can't embed Avenir and dies with "invalid font type".
ggsave_spark(file.path(figures_main_folder, "figure_consort_flowchart.png"),
             fp, width = 9, height = 7.5)
ggsave(file.path(figures_main_folder, "figure_consort_flowchart.pdf"),
       fp, width = 9, height = 7.5, bg = "white", device = grDevices::cairo_pdf)

cat(sprintf("cascade: %d events -> %d cell-days -> minus %d flagged -> %d cases, %s referents\n",
            n_raw_events, n_raw_celldays, n_flagged, n_cases, fmt(n_refs)))
