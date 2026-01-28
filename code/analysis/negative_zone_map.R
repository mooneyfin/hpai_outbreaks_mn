
# ---- Package setup ----
packages <- c(
  "MatchIt", "readr", "sf", "survival", "dlnm", "gtsummary", "kableExtra",
  "tidyr", "tibble", "mgcv", "ggplot2", "car", "patchwork", "corrplot",
  "purrr", "ggspatial", "gt", "flextable", "officer", "conflicted",
  "dplyr", "data.table", "here", "stringr"
)

installed <- packages %in% rownames(installed.packages())
if (any(!installed)) install.packages(packages[!installed])

invisible(lapply(packages, library, character.only = TRUE))

# Resolve conflicts
conflict_prefer("select", "dplyr")
conflict_prefer("filter", "dplyr")

# Project folders 
base_folder <- here()
data_folder <- here("data")
env_folder <- here("data", "EarthEngine_MN")
objects_folder <-here("data", "objects")
figures_folder <- here("output", "figures")
tables_folder <- here("output", "tables")

# Load timeseries clean

timeseries_clean <- readRDS(file.path(objects_folder, "timeseries_clean.rds"))

# Count outbreak zones

outbreak <- timeseries_clean %>%
  group_by(zone_id) %>%
  filter(outbreak_binary == 1)

nrow(no_outbreak_zones)  # number of zones with outbreak

# Count negative control zones

negative_control_zones <- timeseries_clean %>%
  group_by(zone_id, lat, lon, farm_possible, land_cover_label) %>%
  filter(all(outbreak_binary == 0)) %>%
  distinct(zone_id)

nrow(negative_control_zones)  # number of zones with no outbreak

# Read in control zones map (has all of the zones except from those with a case)

control_zones_map <- st_read(file.path(data_folder,"Control_Zones.gpkg"))

# Remove temporal
control_zones_map <- control_zones_map %>%
  group_by(zone_id, geom) %>%
  distinct(zone_id)

nrow(control_zones_map)  # number of zones with no outbreak, including zones without a farm

# Plot quickly
control_zones_sf_proj <- st_transform(control_zones_map, 3857)

ggplot() +
  annotation_map_tile(type = "cartolight", zoomin = 0) +
  # Control zones: outlines
  geom_sf(data = control_zones_sf_proj, fill = NA, color = "lightblue", linewidth = 0.05) 

# What zones have a feedlot but are not cases?

unmatched <- control_zones_map %>%
  anti_join(negative_control_zones, by = "zone_id")

nrow(unmatched)

# Plot all control zones and then negative control zones

negative_control_sf <- control_zones_sf_proj %>%
  filter(zone_id %in% negative_control_zones$zone_id)

ggplot() +
  geom_sf(data = control_zones_sf_proj, fill = "grey90", color = "grey50") +
  geom_sf(data = negative_control_sf, fill = "steelblue", color = "darkblue") +
  theme_minimal()

# EXPORT NEGATIVE CONTROL ZONES FOR MAPPING IN PYTHON

st_write(negative_control_sf, file.path(data_folder, "Negative_Control_Zones.gpkg"))
 


  