# NOTE: This script requires confidential HPAI spillover point-location data
# (not publicly available). The processed grid-cell level output is provided
# in data/ and downstream analyses can be run without this file.

# Validates that outbreak cases in our dataset are truly independent spillover events!
# Cross-references the raw confidential data against zone assignments to check
# for repeat detections at the same farm vs. distinct introductions...

rm(list=ls())
project.folder = paste0(print(here::here()),'/')
source(paste0(project.folder,'create_folder_structure.R'))
source(paste0(functions.folder,'script_initiate.R'))

# ---- Load data ----
spillover_file <- file.path(spillover_data, "hpai-outbreaks-clean-placeholder-032124.csv")
if (!file.exists(spillover_file)) {
  stop("This script needs the confidential spillover data (hpai-outbreaks-clean-placeholder-032124.csv). See 01_data/1b_outcome_data/README.md for how to get it.")
}
hpai_outbreaks <- read_csv(spillover_file)
case_zones <- st_read(file.path(processed_data, "Case_Zones.gpkg"))

# ---- Filter to Minnesota 2022 ----
hpai_mn_2022 <- hpai_outbreaks %>%
  filter(
    state == "MN",
    date_confirmed >= as.Date("2022-01-01"),
    date_confirmed <= as.Date("2022-12-31")
  )

cat("Total MN outbreaks in 2022:", nrow(hpai_mn_2022), "\n")

# ---- Convert to spatial and join to zones ----
hpai_mn_2022_sf <- hpai_mn_2022 %>%
  st_as_sf(coords = c("lon", "lat"), crs = 4326) %>%
  st_transform(st_crs(case_zones))

hpai_with_zones <- hpai_mn_2022_sf %>%
  st_join(case_zones, join = st_intersects)

# Validation: check spatial join
cat("Total outbreaks:", nrow(hpai_mn_2022_sf), "\n")
cat("Outbreaks matched to a zone:", sum(!is.na(hpai_with_zones$zone_id)), "\n")
cat("Outbreaks outside case zones:", sum(is.na(hpai_with_zones$zone_id)), "\n")
## We have left over, because we attribute the outbreak to the feedlot rather than the farm

# ---- Identify unique farm locations within zones ----
hpai_with_zones_flagged <- hpai_with_zones %>%
  filter(!is.na(zone_id)) %>%
  group_by(zone_id) %>%
  mutate(
    n_unique_locations = n_distinct(geometry),
    location_id = as.numeric(factor(as.character(geometry)))
  ) %>%
  group_by(zone_id, location_id) %>%
  mutate(detection_num_at_farm = row_number(date_confirmed)) %>%
  ungroup()

# ---- Define independent spillovers ----
# Keep first detection per farm, OR subsequent detections 8+ weeks later
min_gap_weeks <- 12

hpai_independent <- hpai_with_zones_flagged %>%
  arrange(zone_id, location_id, date_confirmed) %>%
  group_by(zone_id, location_id) %>%
  mutate(
    weeks_since_last = as.numeric(difftime(date_confirmed, lag(date_confirmed), units = "weeks")),
    keep = detection_num_at_farm == 1 | (!is.na(weeks_since_last) & weeks_since_last >= min_gap_weeks)
  ) %>%
  filter(keep) %>%
  ungroup() %>%
  select(-weeks_since_last, -keep)

cat("Independent spillovers:", nrow(hpai_independent), "\n")

# ---- Create case dataset ----
hpai_independent <- hpai_independent %>%
  mutate(case_week = as.numeric(floor(difftime(date_confirmed, as.Date("2022-01-01"), units = "weeks")) + 1))

cases <- hpai_independent %>%
  st_drop_geometry() %>%
  distinct(zone_id, case_week)

cat("Independent cases:", nrow(cases), "\n")

# ---- Validate control-case overlap ----
case_zone_weeks <- cases %>% 
  select(zone_id, case_week)

proposed_controls <- cases %>%
  crossing(lag = 1:4) %>%
  mutate(control_week = case_week - lag) %>%
  filter(control_week >= 1)

overlap <- proposed_controls %>%
  inner_join(case_zone_weeks, by = c("zone_id", "control_week" = "case_week"))

cat("Overlapping controls (retained per study design):", nrow(overlap), "\n")

# ---- Find distribution of cases ----

# First and last case by week
cat("\n=== Case Distribution Summary ===\n")
cat("First case week:", min(hpai_independent$case_week), "\n")
cat("Last case week:", max(hpai_independent$case_week), "\n")

# First and last case by date
cat("\nFirst case date:", as.character(min(hpai_independent$date_confirmed)), "\n")
cat("Last case date:", as.character(max(hpai_independent$date_confirmed)), "\n")

# Total span
date_range <- as.numeric(difftime(max(hpai_independent$date_confirmed), 
                                  min(hpai_independent$date_confirmed), 
                                  units = "days"))
cat("Total span (days):", date_range, "\n")

# Distribution by week
cat("\n=== Cases per Week ===\n")
cases_by_week <- hpai_independent %>%
  st_drop_geometry() %>%
  count(case_week, name = "n_cases") %>%
  arrange(case_week)
print(cases_by_week, n = Inf)

# ---- Examine zones with multiple cases ----
cat("\n=== Zones with Multiple Cases ===\n")

zones_multiple_cases <- hpai_independent %>%
  st_drop_geometry() %>%
  group_by(zone_id) %>%
  summarise(
    n_cases = n(),
    first_date = min(date_confirmed),
    last_date = max(date_confirmed),
    span_days = as.numeric(difftime(max(date_confirmed), min(date_confirmed), units = "days")),
    weeks_with_cases = paste(sort(unique(case_week)), collapse = ", "),
    .groups = "drop"
  ) %>%
  filter(n_cases > 1) %>%
  arrange(desc(n_cases))

cat("Number of zones with >1 case:", nrow(zones_multiple_cases), "\n")
cat("Total cases in multi-case zones:", sum(zones_multiple_cases$n_cases), "\n\n")

print(zones_multiple_cases, n = Inf)

# ---- Summary of temporal clustering in multi-case zones ----
cat("\n=== Temporal Clustering Summary ===\n")

if(nrow(zones_multiple_cases) > 0) {
  cat("Mean span (days) between first and last case in multi-case zones:", 
      round(mean(zones_multiple_cases$span_days), 1), "\n")
  cat("Median span (days):", median(zones_multiple_cases$span_days), "\n")
  cat("Range of spans:", min(zones_multiple_cases$span_days), "-", 
      max(zones_multiple_cases$span_days), "days\n")
  
  # How many had cases on the same day vs spread apart?
  same_day <- sum(zones_multiple_cases$span_days == 0)
  spread_apart <- sum(zones_multiple_cases$span_days > 0)
  cat("\nZones with all cases on same day:", same_day, "\n")
  cat("Zones with cases spread across multiple days:", spread_apart, "\n")
}

# Define season boundaries for 2022
season_lines <- data.frame(
  date = as.Date(c("2022-03-20", "2022-06-21", "2022-09-22", "2022-12-21")),
  label = c("Spring", "Summer", "Fall", "Winter")
)

ggplot(hpai_independent %>% st_drop_geometry(), 
       aes(x = date_confirmed)) +
  # Add histogram
  geom_histogram(binwidth = 7, fill = "black", color = "white") +
  # Add vertical lines at season boundaries
  geom_vline(data = season_lines, aes(xintercept = date), 
             linetype = "dashed", color = "gray40") +
  # Add season labels rotated 90 degrees
  geom_text(data = season_lines, 
            aes(x = date, y = Inf, label = label),
            angle = 90, hjust = 1.1, vjust = -0.3, size = 5, color = "gray30") +
  # X-axis
  scale_x_date(
    limits = c(as.Date("2022-01-01"), as.Date("2022-12-31")),
    date_breaks = "1 month",
    date_labels = "%b",
    expand = c(0, 0)
  ) +
  labs(
    x = "Date Confirmed",
    y = "Number of Cases"
  ) +
  theme_minimal() +
  theme(
    panel.grid.major.x = element_line(color = "gray85", linewidth = 0.3),
    panel.grid.major.y = element_line(color = "gray85", linewidth = 0.3),
    panel.grid.minor = element_blank()
  )

n_distinct(case_crossover_df$stratum_id[case_crossover_df$outbreak_binary == 1])