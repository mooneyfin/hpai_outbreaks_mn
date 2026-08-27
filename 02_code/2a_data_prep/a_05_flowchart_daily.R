rm(list = ls())

#0a.Declare root directory, folder location and load essential stuff
project.folder = paste0(print(here::here()), '/')
source(paste0(project.folder, 'create_folder_structure.R'))
source(paste0(functions.folder, 'script_initiate.R'))

#1.Load cascade objects from a_02 and a_03
summ    <- readRDS(paste0(objects_folder, "spillover_summary.RDS"))
cascade <- readRDS(paste0(objects_folder, "case_crossover_cascade.RDS"))

#2.Extract counts for box labels
total_events     <- summ$total
case_zone_days   <- summ$case_zone_days
case_zones       <- summ$case_zones
multi_event_days <- summ$multi_event_zone_days
non_ind_500m     <- summ$non_independent_500m
cases_final      <- cascade$cases_final
controls_final   <- cascade$controls_final
strata_final     <- cascade$strata_final
strata_dropped   <- cascade$strata_dropped_no_ctrl

#3.Build box content
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

#4.Arrows between consecutive boxes
arrows <- data.frame(
  x    = 0.5,
  xend = 0.5,
  y    = boxes$y[-nrow(boxes)] - 0.04,
  yend = boxes$y[-1]            + 0.04
)

#5.Side annotations
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

#6.Draw flow chart
fp <- ggplot() +
  geom_tile(data = boxes,
            aes(x = x, y = y),
            width = 0.55, height = 0.10,
            fill = "white", color = "grey25", linewidth = 0.4) +
  geom_text(data = boxes,
            aes(x = x, y = y, label = label),
            family = "Avenir", size = 3.4, lineheight = 1.05) +
  geom_segment(data = arrows,
               aes(x = x, xend = xend, y = y, yend = yend),
               arrow = arrow(length = unit(0.18, "cm"), type = "closed"),
               color = "grey25", linewidth = 0.4) +
  geom_text(data = sides,
            aes(x = x, y = y, label = txt),
            family = "Avenir", size = 2.9, hjust = 0,
            color = "grey30", lineheight = 1.05) +
  geom_segment(aes(x = 0.80, xend = 0.80, y = 0.05, yend = 0.98),
               color = "grey80", linewidth = 0.3, linetype = "dotted") +
  scale_x_continuous(limits = c(0, 1), expand = c(0, 0)) +
  scale_y_continuous(limits = c(0, 1), expand = c(0, 0)) +
  labs(
    title    = "Sample-size cascade: time-stratified case-crossover design",
    subtitle = "Highly pathogenic avian influenza (H5N1) spillover into Minnesota poultry, 2022"
  ) +
  theme_spark_map() +
  theme(plot.title.position = "plot",
        plot.margin = margin(15, 15, 15, 15))

#7.Save
ggsave(paste0(figures_main_folder, "figure_daily_flowchart.png"),
       fp, width = 9, height = 11, units = "in", dpi = 300, bg = "white")
ggsave(paste0(figures_main_folder, "figure_daily_flowchart.pdf"),
       fp, width = 9, height = 11, units = "in", bg = "white")
