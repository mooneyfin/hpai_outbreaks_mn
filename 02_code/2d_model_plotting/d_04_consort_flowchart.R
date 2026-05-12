# d_04_consort_flowchart.R
# CONSORT-style sample-size cascade for the daily time-stratified case-crossover
# study of HPAI spillover into MN poultry farms, 2022.
#
# Why this script exists: GeoHealth reviewers asked us (May 2026 revision) to
# make exclusion criteria, sample size, and modeling size explicit. The cascade
# below mirrors the figure they will see in the resubmission.
#
# Inputs : 03_output/3b_model_output/dataframes/spillover_summary.RDS
#          03_output/3b_model_output/dataframes/case_crossover_cascade.RDS
# Outputs: 05_figures/figure_S2_consort_flowchart.{png,pdf}

# 0a. Declare root directory, folder locations, and load packages + helpers
rm(list = ls())
project.folder = paste0(print(here::here()), '/')
source(paste0(project.folder, 'create_folder_structure.R'))
source(paste0(functions.folder, 'script_initiate.R'))

# 0b. Load cascade objects (produced by b_create_timeseries_df.Rmd and
#     c_create_casscrossover_df.Rmd). Each script writes a self-contained list,
#     so this figure regenerates without touching the confidential raw data.
summ    <- readRDS(file.path(objects_folder, "spillover_summary.RDS"))
cascade <- readRDS(file.path(objects_folder, "case_crossover_cascade.RDS"))

# 1a. Pull the key counts used in the box labels — derived from the cascade
#     objects so updates to the upstream pipeline propagate automatically.
total_events       <- summ$total
case_zone_days     <- summ$case_zone_days
case_zones         <- summ$case_zones
multi_event_days   <- summ$multi_event_zone_days
non_ind_500m       <- summ$non_independent_500m
cases_final        <- cascade$cases_final
controls_final     <- cascade$controls_final
strata_final       <- cascade$strata_final
strata_dropped     <- cascade$strata_dropped_no_ctrl

# 2a. Box content (vertical cascade — y positions decrease top to bottom)
boxes <- tibble::tribble(
  ~step, ~x, ~y, ~label,
  1, 0.5, 0.95, sprintf(
    "USDA APHIS confirmed HPAI poultry detections\nMinnesota, 2022\nn = %d events",
    total_events),
  2, 0.5, 0.80, sprintf(
    "Spatial aggregation to 10 km × 1-day grid\n%d events → %d case zone-days\n(events at the same zone-day collapsed)",
    total_events, case_zone_days),
  3, 0.5, 0.63, sprintf(
    "Spatio-temporal independence screen\n< 500 m AND < 28 d from a prior detection\n%d non-independent later events flagged",
    non_ind_500m),
  4, 0.5, 0.46, sprintf(
    "Stratum construction\nzone × year × calendar month × day-of-week\n(Lu et al. 2020 time-stratified case-crossover)"),
  5, 0.5, 0.29, sprintf(
    "Strata dropped (no valid control days in month × DOW slot)\nn = %d", strata_dropped),
  6, 0.5, 0.12, sprintf(
    "Analytic sample\n%d case zone-days  |  %d control zone-days\n%d strata across %d unique zones",
    cases_final, controls_final, strata_final, case_zones)
)

# 2b. Arrows between consecutive boxes
arrows <- data.frame(
  x    = 0.5,
  xend = 0.5,
  y    = boxes$y[-nrow(boxes)] - 0.04,
  yend = boxes$y[-1]            + 0.04
)

# 2c. Side annotations describing what happens at each step (right column)
sides <- data.frame(
  x = 0.82, y = c(0.875, 0.715, 0.545, 0.375, 0.205),
  txt = c(
    "Surveillance events with date_confirmed in 2022",
    sprintf("%d zone-days had >1 event\non the same date", multi_event_days),
    "Strata where the sole event is\nthe non-independent one are removed;\nmulti-event zone-days are retained",
    "Controls = other days in the same\nzone-month-DOW slot that\nare NOT case days",
    "Final controls per case-day:\n3 (months with 4 same-DOWs)\n4 (months with 5 same-DOWs)"
  )
)

# 3a. Draw the flow chart with SPARK Lab figure conventions
#     (theme_spark_map() is the closest existing theme — gives us Avenir, no
#      panel chrome, and Avenir Heavy title/subtitle). The cascade is a
#      non-cartographic ggplot but the typography requirements are identical.
fp <- ggplot2::ggplot() +
  # Box rectangles
  ggplot2::geom_tile(data = boxes,
                     ggplot2::aes(x = x, y = y),
                     width = 0.55, height = 0.10,
                     fill = "white", color = "grey25", linewidth = 0.4) +
  # Box labels
  ggplot2::geom_text(data = boxes,
                     ggplot2::aes(x = x, y = y, label = label),
                     family = "Avenir", size = 3.4, lineheight = 1.05) +
  # Cascade arrows
  ggplot2::geom_segment(data = arrows,
                        ggplot2::aes(x = x, xend = xend, y = y, yend = yend),
                        arrow = ggplot2::arrow(length = ggplot2::unit(0.18, "cm"), type = "closed"),
                        color = "grey25", linewidth = 0.4) +
  # Right-column side annotations
  ggplot2::geom_text(data = sides,
                     ggplot2::aes(x = x, y = y, label = txt),
                     family = "Avenir", size = 2.9, hjust = 0,
                     color = "grey30", lineheight = 1.05) +
  # Right-of-cascade vertical guide
  ggplot2::geom_segment(ggplot2::aes(x = 0.80, xend = 0.80, y = 0.05, yend = 0.98),
                        color = "grey80", linewidth = 0.3, linetype = "dotted") +
  ggplot2::scale_x_continuous(limits = c(0, 1), expand = c(0, 0)) +
  ggplot2::scale_y_continuous(limits = c(0, 1), expand = c(0, 0)) +
  ggplot2::labs(
    title    = "Sample-size cascade: time-stratified case-crossover design",
    subtitle = "Highly pathogenic avian influenza (H5N1) spillover into Minnesota poultry, 2022"
  ) +
  theme_spark_map() +
  ggplot2::theme(plot.title.position = "plot",
                 plot.margin = ggplot2::margin(15, 15, 15, 15))

# 4a. Export (300 DPI; PNG for the manuscript figures folder, PDF for journals
#     that prefer vector). Width/height chosen for a portrait single-column page.
ggplot2::ggsave(file.path(figures_folder, "figure_S2_consort_flowchart.png"),
                fp, width = 9, height = 11, units = "in", dpi = 300, bg = "white")
ggplot2::ggsave(file.path(figures_folder, "figure_S2_consort_flowchart.pdf"),
                fp, width = 9, height = 11, units = "in", bg = "white")

cat("Wrote CONSORT flowchart to ",
    file.path(figures_folder, "figure_S2_consort_flowchart.png"), "\n", sep = "")
