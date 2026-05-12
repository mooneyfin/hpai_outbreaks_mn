rm(list = ls())

# 0a. Load packages
library(here)

# 0b. Declare directories — SPARK Lab NYC layout (May 2026 revision)
project.folder = paste0(print(here::here()), '/')

# Code (02_code/)
code.folder              = paste0(project.folder, "02_code/")
  functions.folder       = paste0(code.folder,    "20_functions/")
  packages.folder        = paste0(code.folder,    "20_functions/")
  file.locations.folder  = paste0(code.folder,    "20_functions/")
  data.prep.folder       = paste0(code.folder,    "2a_data_prep/")
  data.exploration.folder= paste0(code.folder,    "2b_data_exploration/")
  models.code.folder     = paste0(code.folder,    "2c_models/")
  model.processing.folder= paste0(code.folder,    "2d_model_plotting/")

# Data (01_data/)
data.folder              = paste0(project.folder, "01_data/")
  exposure.data.folder   = paste0(data.folder,    "1a_exposure_data/")
  outcome.data.folder    = paste0(data.folder,    "1b_outcome_data/")
  supportive.data.folder = paste0(data.folder,    "1c_supportive_datasets/")
    env_data             = paste0(exposure.data.folder,   "meteorological_data/")
    bird_abundance_data  = paste0(exposure.data.folder,   "bird_abundance_data/")
    spillover_data       = outcome.data.folder
    feedlot_data         = paste0(supportive.data.folder, "feedlot_data/")
    spatial_data         = paste0(supportive.data.folder, "spatial_data/")

# Output (03_output/)
output.folder            = paste0(project.folder, "03_output/")
  eda.output.folder      = paste0(output.folder,  "3a_eda_output/")
  model.output.folder    = paste0(output.folder,  "3b_model_output/")
    objects_folder       = paste0(model.output.folder, "dataframes/")  # analysis-ready RDS
    models_folder        = paste0(model.output.folder, "models/")      # fitted models
    processed_data       = eda.output.folder                            # diagnostic outputs (gpkg/csv)

# Tables, figures, drafts
tables_folder            = paste0(project.folder, "04_tables/")
figures_folder           = paste0(project.folder, "05_figures/")
literature.folder        = paste0(project.folder, "06_literature/")
drafts.folder            = paste0(project.folder, "07_drafts/")

# 1a. Create any missing folders
all_folders = c(
  code.folder, functions.folder, data.prep.folder, data.exploration.folder,
  models.code.folder, model.processing.folder,
  data.folder, exposure.data.folder, outcome.data.folder, supportive.data.folder,
  env_data, bird_abundance_data, feedlot_data, spatial_data,
  output.folder, eda.output.folder, model.output.folder,
  objects_folder, models_folder,
  tables_folder, figures_folder, literature.folder, drafts.folder
)

for (f in all_folders) {
  if (!dir.exists(f)) dir.create(f, recursive = TRUE)
}
