# Run the full multi-state analysis end to end: data prep (a_04 -> a_01 -> a_02 -> a_03)
# then the gridded case-crossover models (c_05). Each Rmd is executed in its OWN Rscript
# subprocess, because each script begins with rm(list = ls()) — sourcing them in a single
# session would wipe this runner's own loop state. Run from the project root:
#   Rscript 02_code/run_multistate_pipeline.R

files <- c(
  "02_code/2a_data_prep/a_04_zone_bird_abundance_multistate.Rmd",
  "02_code/2a_data_prep/a_01_join_data_multistate.Rmd",
  "02_code/2a_data_prep/a_02_create_timeseries_df_multistate.Rmd",
  "02_code/2a_data_prep/a_03_create_casecrossover_df_multistate.Rmd",
  "02_code/2c_models/c_05_casecrossover_multistate_daily.Rmd",
  "02_code/2c_models/c_06_casecrossover_multistate_weekly.Rmd"
)

run_one <- function(f) {
  code <- sprintf(
    "suppressWarnings(suppressMessages({r<-knitr::purl('%s',output=tempfile(fileext='.R'),quiet=TRUE);source(r,echo=FALSE)}))",
    f)
  system2("Rscript", args = c("-e", shQuote(code)), stdout = "", stderr = "")
}

for (f in files) {
  cat("\n##########", f, "##########\n")
  t0 <- proc.time()[["elapsed"]]
  status <- run_one(f)
  if (!is.null(status) && status != 0) {
    cat("!!! FAILED (exit", status, "):", f, "\n"); quit(status = 1)
  }
  cat(sprintf(">>> done in %.1f s\n", proc.time()[["elapsed"]] - t0))
}
cat("\n########## FULL ANALYSIS COMPLETE ##########\n")
