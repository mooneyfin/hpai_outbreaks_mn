#1.Centralized file path declarations
source(paste0(file.locations.folder, 'file_locations.R'))

#2.Load packages used across all analyses
source(paste0(packages.folder, 'packages_to_load.R'))

#3.Load shared helper functions (figure theme, palettes, etc.)
source(paste0(functions.folder, 'functions.R'))

#4. Register SPARK Lab font variants so ragg/systemfonts can find them
library(systemfonts)
if (requireNamespace("systemfonts", quietly = TRUE)) {
  variants_needed <- list(
    list(name = "Avenir Heavy",  weight = "heavy"),
    list(name = "Avenir Medium", weight = "medium"),
    list(name = "Avenir Book",   weight = "normal")
  )
  for (v in variants_needed) {
    try(systemfonts::register_variant(
      name   = v$name,
      family = "Avenir",
      weight = v$weight
    ), silent = TRUE)
  }
}