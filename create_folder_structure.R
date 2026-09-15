rm(list = ls())

#0.Load packages
library(here)

#1a.Declare directories (can add to over time)
project.folder <- paste0(print(here::here()), '/')

data.folder <- paste0(project.folder, "01_data/")
  exposure.data.folder   <- paste0(data.folder, "1a_exposure_data/")
  outcome.data.folder    <- paste0(data.folder, "1b_outcome_data/")
  supportive.data.folder <- paste0(data.folder, "1c_supportive_datasets/")
    env_data             <- paste0(exposure.data.folder,   "meteorological_data/")
    bird_abundance_data  <- paste0(exposure.data.folder,   "bird_abundance_data/")
    spillover_data       <- outcome.data.folder
    feedlot_data         <- paste0(supportive.data.folder, "feedlot_data/")
    spatial_data         <- paste0(supportive.data.folder, "spatial_data/")

code.folder <- paste0(project.folder, "02_code/")
  data.prep.folder        <- paste0(code.folder, "2a_data_prep/")
  data.exploration.folder <- paste0(code.folder, "2b_data_exploration/")
  models.code.folder      <- paste0(code.folder, "2c_models/")
  model.processing.folder <- paste0(code.folder, "2d_model_plotting/")
  functions.folder        <- paste0(code.folder, "20_functions/")
  packages.folder         <- functions.folder
  file.locations.folder   <- functions.folder

output.folder <- paste0(project.folder, "03_output/")
  eda.output.folder   <- paste0(output.folder, "3a_eda_output/")
  model.output.folder <- paste0(output.folder, "3b_model_output/")
    objects_folder   <- paste0(model.output.folder, "dataframes/")
    models_folder    <- paste0(model.output.folder, "models/")
    processed_data   <- eda.output.folder

tables_folder  <- paste0(project.folder, "04_tables/")
  tables_main_folder         <- paste0(tables_folder,  "main/")          # Table 1
  tables_supplement_folder   <- paste0(tables_folder,  "supplement/")    # Tables S2-S9
  tables_specs_folder        <- paste0(tables_folder,  "table_specs/")   # display specs the brief renders from
  tables_archive_folder      <- paste0(tables_folder,  "archive/")
figures_folder <- paste0(project.folder, "05_figures/")
  figures_main_folder         <- paste0(figures_folder, "main/")          # Figures 1-4
  figures_supplement_folder   <- paste0(figures_folder, "supplement/")    # Figures S1-S4
  figures_intermediate_folder <- paste0(figures_folder, "intermediate/")  # map and curve that compose Figure 1
  figures_archive_folder      <- paste0(figures_folder, "archive/")
lit.folder     <- paste0(project.folder, "06_literature/")
drafts.folder  <- paste0(project.folder, "07_drafts/")

#1b.Identify list of folder locations which have just been created above
folders.names <- grep("\\.folder$|_folder$|_data$", names(.GlobalEnv), value = TRUE)

#1c.Create function to create folders
create_folders <- function(name) {
  ifelse(!dir.exists(get(name)), dir.create(get(name), recursive = TRUE), FALSE)
}

#1d.Create the folders named above
invisible(lapply(folders.names, create_folders))
