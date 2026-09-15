## ============================================================
## Sample-size cascade — weekly time-stratified case-crossover
## (HPAI / H5N1 poultry spillover in Minnesota, 2022)
## ============================================================

rm(list = ls())

# 0a. Declare root directory, folder location and load essential stuff
project.folder <- paste0(print(here::here()), "/")
source(paste0(project.folder, "create_folder_structure.R"))
source(paste0(functions.folder, "script_initiate.R"))

suppressPackageStartupMessages({
  library(ggplot2); library(sf); library(data.table); library(tibble); library(grid)
})

# 1. Load source data, fishnet, and analytic bundle
spillover_csv <- fread(paste0(spillover_data,
                              "hpai-outbreaks-clean-placeholder-032124.csv"))
spillover_csv[, year := as.integer(substr(date_confirmed, 1, 4))]
mn22 <- spillover_csv[state == "MN" & year == 2022]

fishnet <- st_read(paste0(spatial_data, "fishnet_10km.shp"), quiet = TRUE)
summ    <- readRDS(paste0(objects_folder, "spillover_summary.RDS"))
ab      <- readRDS(paste0(objects_folder, "weekly_primary_ab.RDS"))
cc_A    <- ab$cc_A; setDT(cc_A)

# Spatial join — events to fishnet zones (10 km grid)
events_sf <- st_as_sf(mn22, coords = c("lon", "lat"), crs = 4326) |>
  st_transform(st_crs(fishnet))
joined <- st_join(events_sf, fishnet["zone_id_ne"], join = st_within) |>
  st_drop_geometry()
setDT(joined)

# 2. Counts (computed live so the cascade reconciles end-to-end)
total_events       <- nrow(mn22)                                        # 110
total_farms        <- uniqueN(mn22$premises)                            # 110
total_commercial   <- uniqueN(mn22[herd_flock_type == "Commercial", premises])  # 81
total_backyard     <- uniqueN(mn22[herd_flock_type == "Backyard",   premises])  # 29
total_zones        <- uniqueN(joined$zone_id_ne)                        # 76
cslt_removed       <- mn22[link == "CSLT",  .N]                         # 20
noseq_removed      <- mn22[link == "NOSEQ", .N]                         # 2
events_retained    <- mn22[link %in% c("IND", "INR"), .N]               # 88
farms_retained     <- mn22[link %in% c("IND", "INR"), uniqueN(premises)]
zones_retained     <- uniqueN(joined[link %in% c("IND", "INR"), zone_id_ne])  # 70
non_independent    <- summ$non_independent_500m                         # 2
cases_final        <- sum(cc_A$outbreak_binary == 1)                    # 85
controls_final     <- sum(cc_A$outbreak_binary == 0)                    # 340
strata_final       <- uniqueN(cc_A$stratum_id)                          # 85
zones_final        <- uniqueN(cc_A$zone_id)                             # 70
multi_per_week     <- events_retained - cases_final                     # 3

# 3. Boxes — 5 sections, body text only (headers rendered separately)
boxes <- tribble(
  ~step, ~x, ~y,    ~fill,     ~header,             ~label,
  1, 0.5, 0.880, "#F4F4F4", "Data source", sprintf(
    "USDA APHIS confirmed HPAI poultry detections\nMinnesota, 2022\n%d events at %d unique commercial (%d) and backyard (%d) farms\nacross %d unique 10 km zones",
    total_events, total_farms, total_commercial, total_backyard, total_zones),
  2, 0.5, 0.685, "#FBEFE4", "Exclusions", sprintf(
    "%d common-source lateral transfer (CSLT) events removed\n%d events without sequence data (NOSEQ) removed\n%d non-independent later events (<500 meters, <28 days from prior in same zone)",
    cslt_removed, noseq_removed, non_independent),
  3, 0.5, 0.490, "#E8F0F5", "Spatial aggregation", sprintf(
    "%d retained events at %d farms collapsed onto %d unique 10 km zones",
    events_retained, farms_retained, zones_retained),
  4, 0.5, 0.295, "#E8F0F5", "Conditional logistic regression", sprintf(
    "outbreak_binary = 1 if any case in zone-week (ISO week, Monday start)\n%d unique case-weeks (%d zones with >1 event in same week)\n1:4 case-control matching: case-week W with controls at W-1, W-2, W-3, W-4 (same zone)",
    cases_final, multi_per_week),
  5, 0.5, 0.100, "#E1EBD8", "Analytic sample", sprintf(
    "%d case-weeks  |  %d control-weeks\n%d strata across %d unique zones",
    cases_final, controls_final, strata_final, zones_final)
)

# Header sits just above each box; box body text stays centered inside the box
box_half_height <- 0.062
header_offset   <- 0.020                           # distance above box top
boxes$header_y  <- boxes$y + box_half_height + header_offset

# 4. Arrows — connect between sections; head terminates above next header
arrow_top_gap    <- 0.010
arrow_bottom_gap <- 0.012
arrows <- data.frame(
  x    = 0.5, xend = 0.5,
  y    = boxes$y[-nrow(boxes)] - box_half_height - arrow_top_gap,
  yend = boxes$header_y[-1]      + arrow_bottom_gap
)

# 5. Draw
fp <- ggplot() +
  # boxes
  geom_tile(data = boxes,
            aes(x = x, y = y, fill = fill),
            width = 0.88, height = box_half_height * 2,
            color = "grey25", linewidth = 0.6) +
  scale_fill_identity() +
  # arrows (rendered first so header text can sit cleanly above the arrowhead)
  geom_segment(data = arrows,
               aes(x = x, xend = xend, y = y, yend = yend),
               arrow = arrow(angle = 22, length = unit(0.36, "cm"),
                             type = "closed"),
               color = "grey20", linewidth = 0.9, lineend = "round") +
  # bold section headers, positioned just above each box
  geom_text(data = boxes,
            aes(x = x, y = header_y, label = header),
            family = "Avenir", fontface = "bold", size = 4.6,
            color = "grey10") +
  # box body text — centred inside each box
  geom_text(data = boxes,
            aes(x = x, y = y, label = label),
            family = "Avenir", size = 4.0, lineheight = 1.20,
            color = "grey15") +
  scale_x_continuous(limits = c(0, 1), expand = c(0, 0)) +
  scale_y_continuous(limits = c(0, 1), expand = c(0, 0)) +
  labs(
    title = "Sample-size: Weekly time-stratified case-crossover"
  ) +
  theme_void(base_family = "Avenir") +
  theme(
    plot.title          = element_text(family = "Avenir", face = "bold",
                                        size = 15, color = "grey10",
                                        hjust = 0.5,
                                        margin = margin(b = 24)),
    plot.title.position = "plot",
    plot.margin         = margin(24, 24, 24, 24)
  )

# 6. Save
ggsave(paste0(figures_main_folder, "figure_weekly_flowchart.png"),
       fp, width = 10, height = 13, units = "in", dpi = 300, bg = "white")

cat("Saved: figure_weekly_flowchart.png  (",
    cases_final, " cases / ", controls_final, " controls / ",
    strata_final, " strata across ", zones_final, " zones; source: ",
    total_events, " events at ", total_farms, " farms in ", total_zones,
    " zones)\n", sep = "")
