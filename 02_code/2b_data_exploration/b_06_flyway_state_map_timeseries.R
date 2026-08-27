## ============================================================
## EDA — HPAI commercial poultry detections, 7-state northern
## Mississippi flyway (MN, WI, MI, IA, IL, IN, OH), 2022
## Output: 05_figures/figure_flyway_state_panel.png
## CSLT (common-source lateral transfer) and NOSEQ events excluded.
## ============================================================

rm(list = ls())

# 0a. Declare root directory, folder location and load essential stuff
project.folder <- paste0(print(here::here()), "/")
source(paste0(project.folder, "create_folder_structure.R"))
source(paste0(functions.folder, "script_initiate.R"))

suppressPackageStartupMessages({
  library(ggplot2); library(sf); library(data.table); library(lubridate)
  library(viridisLite); library(patchwork)
})

# 1. Load outbreak data, filter to 7-state, 2022, post-CSLT panel
d <- fread(paste0(spillover_data, "hpai-outbreaks-clean-placeholder-032124.csv"))
d[, year := as.integer(substr(date_confirmed, 1, 4))]

flyway_states <- c("MN", "WI", "MI", "IA", "IL", "IN", "OH")
flyway_full   <- c("Minnesota", "Wisconsin", "Michigan", "Iowa",
                   "Illinois",  "Indiana",   "Ohio")

events <- d[year       == 2022 &
            state      %in% flyway_states &
            link       %in% c("IND", "INR")]
events[, date := as.Date(date_confirmed)]

# State-level counts, ordered DESCENDING so the highest count anchors the
# dark end of the viridis ramp (Minnesota → dark, Illinois → light).
state_counts <- events[, .N, by = state][order(-N)]
state_counts[, state_full := factor(state, levels = flyway_states,
                                    labels = flyway_full)]
state_counts[, state_full := factor(state_full,
                                    levels = as.character(state_full))]

# 2. Build the shared viridis palette: 7 evenly-spaced colours, one per state
# (ordered by descending count), so each state has a visually distinct hue
# rather than crowding the lower-count states into similar greens.
# direction = -1 anchors MN (highest count) at the dark-purple end.
n_states      <- nrow(state_counts)
state_palette <- setNames(
  viridis(n_states, end = 0.92),                     # default direction: dark purple → yellow
  as.character(state_counts$state_full)              # state_counts is ordered descending by count
)

# 3. State polygons from TIGER 2018 (mirrors the GEE multistate script)
states_url <- "https://www2.census.gov/geo/tiger/GENZ2022/shp/cb_2022_us_state_5m.zip"
local_dir  <- paste0(spatial_data, "cb_2022_us_state_5m/")
local_shp  <- paste0(local_dir,  "cb_2022_us_state_5m.shp")
if (!file.exists(local_shp)) {
  dir.create(local_dir, recursive = TRUE, showWarnings = FALSE)
  tmpzip <- tempfile(fileext = ".zip")
  download.file(states_url, tmpzip, mode = "wb", quiet = TRUE)
  unzip(tmpzip, exdir = local_dir)
}
states_sf <- st_read(local_shp, quiet = TRUE)
albers <- "+proj=aea +lat_1=29.5 +lat_2=45.5 +lat_0=37.5 +lon_0=-96 +datum=NAD83"

flyway_sf <- states_sf[states_sf$STUSPS %in% flyway_states, ] |>
  st_transform(albers)
flyway_sf <- merge(flyway_sf, state_counts[, .(state, state_full, N)],
                   by.x = "STUSPS", by.y = "state", all.x = TRUE)
flyway_sf$state_full <- factor(flyway_sf$state_full,
                                levels = levels(state_counts$state_full))

flyway_sf$centroid <- st_centroid(st_geometry(flyway_sf))
centroids <- st_coordinates(flyway_sf$centroid)
flyway_sf$lon_c <- centroids[, 1]
flyway_sf$lat_c <- centroids[, 2]
bbox <- st_bbox(flyway_sf)

# 4. Map: fill by state (categorical), counts shown as in-state labels
p_map <- ggplot(flyway_sf) +
  geom_sf(aes(fill = state_full), color = "white", linewidth = 0.4) +
  geom_text(aes(x = lon_c, y = lat_c,
                label = paste0(STUSPS, "\nn = ", N)),
            family = "Avenir", fontface = "bold",
            size = 4.0, color = "white", lineheight = 1.05) +
  scale_fill_manual(values = state_palette, guide = "none") +
  coord_sf(xlim = c(bbox["xmin"], bbox["xmax"]),
           ylim = c(bbox["ymin"], bbox["ymax"]),
           expand = TRUE) +
  theme_void(base_family = "Avenir") +
  theme(plot.margin = margin(8, 8, 8, 8))

# 5. Timeseries: weekly stacked bars, same palette
events[, week_start := floor_date(date, unit = "week", week_start = 1)]
events[, state_full := factor(state, levels = flyway_states, labels = flyway_full)]
events[, state_full := factor(state_full,
                              levels = levels(state_counts$state_full))]
weekly <- events[, .N, by = .(week_start, state_full)]

# Astronomical seasons (Northern Hemisphere) for 2022
seasons <- data.frame(
  date  = as.Date(c("2022-03-20", "2022-06-21",
                    "2022-09-22", "2022-12-21")),
  label = c("Spring", "Summer", "Fall", "Winter")
)

p_ts <- ggplot(weekly, aes(x = week_start, y = N, fill = state_full)) +
  geom_vline(data = seasons, aes(xintercept = as.numeric(date)),
             linetype = "dashed", color = "grey55", linewidth = 0.4) +
  geom_col(width = 6, color = "white", linewidth = 0.15) +
  geom_text(data = seasons,
            aes(x = date, y = Inf, label = label),
            inherit.aes = FALSE, family = "Avenir", size = 3.0,
            color = "grey35", vjust = 1.6, hjust = -0.05) +
  scale_fill_manual(values = state_palette,
                     name   = "Spillover events excl. CSLT",
                     drop   = FALSE) +
  scale_x_date(date_breaks = "1 month", date_labels = "%b") +
  scale_y_continuous(breaks = seq(0, 20, by = 2)) +
  labs(x = "Week of confirmed detection (2022)", y = "Weekly count") +
  theme_classic(base_family = "Avenir") +
  theme(
    axis.text          = element_text(size = 10),
    axis.title         = element_text(size = 11),
    legend.position    = "right",
    panel.grid.major.y = element_line(color = "grey92", linewidth = 0.3),
    plot.margin        = margin(8, 8, 8, 8)
  )

# 6. Patchwork: map | timeseries, legend collected on the right
panel <- (p_map | p_ts) +
  plot_layout(widths = c(1, 1.7), guides = "collect") &
  theme(legend.position = "right",
        legend.title    = element_text(family = "Avenir", face = "bold",
                                        size = 10, color = "grey15"),
        legend.text     = element_text(family = "Avenir", size = 9))

ggsave(paste0(figures_eda_folder, "figure_flyway_state_panel.png"),
       panel, width = 16, height = 6.5, units = "in",
       dpi = 300, bg = "white")

cat("Saved: figure_flyway_state_panel.png\n")
