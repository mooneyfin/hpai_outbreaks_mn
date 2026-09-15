############ SUPPLEMENT — WHAT HAPPENED TO SURFACE RUNOFF ############
# Runoff was one of the strongest predictors in the original submission. It is not in the
# revised primary model, and a reader is entitled to know why. This tests it every way we
# reasonably can and reports all of them in one table.
#
# Three models on the current specification: runoff on its own, runoff with the other
# non-water exposures but no competing water term, and runoff added to the primary model.
# Null in all three, on all three panels.

rm(list = ls())
project.folder = paste0(print(here::here()), '/')
source(paste0(project.folder, 'create_folder_structure.R'))
source(paste0(functions.folder, 'script_initiate.R'))
suppressMessages({library(INLA); library(Matrix); library(flextable); library(officer)})
INLA::inla.setOption(num.threads = 1)
source(paste0(functions.folder, 'inla_dlnm_helpers.R'))

PREC <- 4; NS2 <- list(fun = "ns", df = 2); LIN <- list(fun = "lin")
cell <- function(e) sprintf("%.2f (%.2f, %.2f)", e[["OR"]], e[["low"]], e[["high"]])

fit_get <- function(d, mv, target, av = NULL) {
  av <- av %||% c(setNames(rep(list(NS2), length(mv) - 1), head(mv, -1)),
                  list(anseriformes = LIN))
  o <- assemble_cc_dlnm(d, met_vars = mv, weekly_vars = character(0), max_lag = 28,
                        argvar_by_var = av, arglag = list(fun = "ns", df = 2))
  f <- suppressWarnings(fit_cc_inla_dlnm(o, fixed_prec = PREC))
  list(waic = f$waic$waic,
       est  = cell(met_lag_effect(f, o, target, 0:7)),
       bird = cell(met_lag_effect(f, o, "anseriformes", 22:28)))
}

rows <- list()
add <- function(model, panel, r) {
  rows[[length(rows) + 1L]] <<- data.table(
    Model = model, Panel = panel, WAIC = round(r$waic, 1),
    `Runoff 0-7 days` = r$est, `Anseriformes 22-28 days` = r$bird)
}

# Current specification throughout: time-stratified one-month strata, post-only 7-day
# washout, 90-day anomalies, N(0, 0.5^2). Only the water terms change.
D <- setDT(readRDS(paste0(objects_folder, "case_crossover_df_shock_post7.RDS")))
if (anyDuplicated(names(D))) D <- D[, which(!duplicated(names(D))), with = FALSE]

MODELS <- list(
  `Runoff alone`                         = c("runoff_shock"),
  `Runoff, no other water terms`         = c("temperature_shock", "runoff_shock",
                                             "wind_speed_shock"),
  `Runoff added to the primary model`    = c("temperature_shock", "precipitation_shock",
                                             "soil_moisture_shock", "wind_speed_shock",
                                             "runoff_shock"))
for (pn in c("Minnesota", "Other flyway states", "Pooled")) {
  d <- switch(pn, Minnesota = D[state == "Minnesota"],
              `Other flyway states` = D[state != "Minnesota"], copy(D))
  for (nm in names(MODELS))
    add(nm, pn, fit_get(d, c(MODELS[[nm]], "anseriformes"), "runoff_shock"))
}

out <- rbindlist(rows)
saveRDS(out, paste0(objects_folder, "tableS5_runoff.RDS"))
is_sig <- function(x) {
  m <- regmatches(x, regexec("\\(([0-9.]+), ([0-9.]+)\\)", x))
  vapply(m, function(z) length(z) == 3 && (as.numeric(z[2]) > 1 | as.numeric(z[3]) < 1), TRUE)
}
ft <- flextable(out) |>
  add_header_row(values = c("", "Odds ratio (95% credible interval) per +0.5 SD"),
                 colwidths = c(3, 2)) |>
  fontsize(size = 10, part = "all") |> autofit() |>
  align(align = "center", part = "all") |> theme_booktabs()
for (j in c("Runoff 0-7 days", "Anseriformes 22-28 days")) {
  hit <- which(is_sig(out[[j]])); if (length(hit)) ft <- bold(ft, i = hit, j = j)
}
read_docx() |> body_add_flextable(ft) |>
  print(target = file.path(tables_main_folder, "TableS5_runoff.docx"))

print(out)
cat("\nwrote TableS5_runoff.docx\n")
