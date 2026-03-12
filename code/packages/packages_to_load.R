############ Packages for HPAI case-crossover analysis ############

# Core packages used across the pipeline
list.of.packages = c(
  'conflicted', 'data.table', 'dplyr', 'here', 'readr', 'tidyr', 'tibble',
  'purrr', 'stringr', 'forcats',
  'sf', 'ggspatial', 'leaflet',
  'survival', 'dlnm', 'mgcv', 'clogitL1', 'MatchIt',
  'ggplot2', 'patchwork', 'corrplot', 'cowplot', 'viridis',
  'gtsummary', 'flextable', 'officer', 'kableExtra', 'gt', 'table1',
  'car', 'boot', 'broom', 'broom.mixed', 'gratia', 'ggeffects',
  'lme4', 'lmerTest', 'emmeans', 'clubSandwich'
)

# Install anything missing
new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]
if(length(new.packages)) install.packages(new.packages, repos = "https://cloud.r-project.org")

# Load all
lapply(list.of.packages, require, character.only = TRUE)

# Resolve namespace conflicts
conflict_prefer("select", "dplyr")
conflict_prefer("filter", "dplyr")
conflict_prefer("lag", "dplyr")
conflict_prefer("label", "table1")
