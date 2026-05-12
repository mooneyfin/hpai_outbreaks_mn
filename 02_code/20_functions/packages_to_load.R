############ Packages for HPAI case-crossover analysis ############
# SPARK Lab NYC pattern: auto-install missing packages and load all silently.
# Note: rgee is NOT in the core list — it requires a one-time Python setup
# (rgee::ee_install()) and GEE auth. Load it on demand inside the rgee
# acquisition script only.

# 1a. Add new packages here, as necessary
list.of.packages = c(
  'conflicted', 'data.table', 'dplyr', 'here', 'readr', 'tidyr', 'tibble',
  'purrr', 'stringr', 'forcats', 'lubridate',
  'sf', 'ggspatial', 'leaflet',
  'survival', 'dlnm', 'mgcv', 'clogitL1', 'MatchIt', 'gnm',
  'ggplot2', 'patchwork', 'corrplot', 'cowplot', 'viridis',
  'gtsummary', 'flextable', 'officer', 'kableExtra', 'gt', 'table1',
  'car', 'boot', 'broom', 'broom.mixed', 'gratia', 'ggeffects',
  'lme4', 'lmerTest', 'emmeans', 'clubSandwich',
  'e1071'
)

# 1b. Check if list of packages is installed. If not, install missing ones.
new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]
if (length(new.packages)) {
  invisible(install.packages(new.packages, repos = "https://cloud.r-project.org"))
}

# 1c. Load packages
invisible(lapply(list.of.packages, require, character.only = TRUE, quietly = TRUE))

# 1d. Resolve namespace conflicts (data.table::select vs dplyr::select, etc.)
conflict_prefer("select", "dplyr",  quiet = TRUE)
conflict_prefer("filter", "dplyr",  quiet = TRUE)
conflict_prefer("lag",    "dplyr",  quiet = TRUE)
conflict_prefer("label",  "table1", quiet = TRUE)
