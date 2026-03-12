rm(list=ls())

# 0a Load Packages
library(here)

# 1a Declare directories
project.folder = paste0(print(here::here()),'/')

  # Code folders
  code.folder = paste0(project.folder, "code/")
    data.prep.folder = paste0(code.folder, "data_prep/")
    data.exploration.folder = paste0(code.folder, "data_exploration/")
    functions.folder = paste0(code.folder, "functions/")
    packages.folder = paste0(code.folder, "packages/")
    models.code.folder = paste0(code.folder, "models/")
    model.processing.folder = paste0(code.folder, "model_processing/")

  # Data folders
  data.folder = paste0(project.folder, "data/")
    file.locations.folder = paste0(data.folder, "file_locations/")
    objects.folder = paste0(data.folder, "objects/")
    objects_folder = paste0(data.folder, "objects/dataframes/")
    models_folder = paste0(data.folder, "objects/models/")
    processed_data = paste0(data.folder, "processed/")
    spillover_data = paste0(data.folder, "confidential/")
    feedlot_data = paste0(data.folder, "raw/feedlot_data/")
    env_data = paste0(data.folder, "raw/meteorological_data/")
    bird_abundance_data = paste0(data.folder, "raw/bird_abundance_data/")
    spatial_data = paste0(data.folder, "raw/spatial_data/")

  # Output folders
  output.folder = paste0(project.folder, "output/")
  figures_folder = paste0(project.folder, "figures/")
  tables_folder = paste0(project.folder, "tables/")
  reports.folder = paste0(project.folder, "reports/")

# 1b Create all folders if they don't exist
all_folders = c(
  code.folder, data.prep.folder, data.exploration.folder, functions.folder,
  packages.folder, models.code.folder, model.processing.folder,
  data.folder, file.locations.folder, objects.folder, objects_folder,
  models_folder, processed_data, spillover_data, feedlot_data, env_data,
  bird_abundance_data, spatial_data, output.folder, figures_folder,
  tables_folder, reports.folder
)

for (f in all_folders) {
  if (!dir.exists(f)) dir.create(f, recursive = TRUE)
}
