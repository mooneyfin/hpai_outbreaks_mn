#1a.Add new packages here, as necessary
list.of.packages = c(
  'conflicted', 'data.table', 'dplyr', 'here', 'readr', 'tidyr', 'tibble',
  'purrr', 'stringr', 'forcats', 'lubridate',
  'sf', 'ggspatial', 'leaflet', 'rgee',
  'survival', 'dlnm', 'mgcv', 'clogitL1', 'MatchIt', 'gnm',
  'ggplot2', 'patchwork', 'corrplot', 'cowplot', 'viridis', 'RColorBrewer',
  'gtsummary', 'flextable', 'officer', 'kableExtra', 'gt', 'table1',
  'car', 'boot', 'broom', 'broom.mixed', 'gratia', 'ggeffects',
  'lme4', 'lmerTest', 'emmeans', 'clubSandwich',
  'e1071', "geojsonio"
)

#1b.Check if list of packages is installed. If not, install missing ones
new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]
if (length(new.packages)) invisible(install.packages(new.packages, repos = "https://cloud.r-project.org"))

#1c.Load packages
invisible(lapply(list.of.packages, require, character.only = TRUE, quietly = TRUE))

#1c-INLA. INLA is NOT on CRAN, so it is loaded separately and never auto-installed.
#         Used by the Bayesian case-crossover DLNM (c_11 + 20_functions/inla_dlnm_helpers.R).
#         Install once with:
#           install.packages("INLA",
#             repos = c(getOption("repos"),
#                       INLA = "https://inla.r-inla-download.org/R/stable"), dep = TRUE)
if (!requireNamespace("INLA", quietly = TRUE)) {
  message("INLA not installed; Bayesian DLNM (c_11) unavailable. See packages_to_load.R for install command.")
} else {
  suppressMessages(require(INLA, quietly = TRUE))
}

#1d.Resolve namespace conflicts
conflict_prefer("select", "dplyr",      quiet = TRUE)
conflict_prefer("filter", "dplyr",      quiet = TRUE)
conflict_prefer("lag",    "dplyr",      quiet = TRUE)
conflict_prefer("label",  "table1",     quiet = TRUE)
conflict_prefer("first",  "data.table", quiet = TRUE)
conflict_prefer("last",   "data.table", quiet = TRUE)
