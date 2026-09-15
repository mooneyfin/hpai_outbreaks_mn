############ INTERACTION VARIANTS — the sensitivities the Figure 4 text cites ############
# c_07 is the single pre-stated interaction model. Its robustness is stated in the manuscript
# in one clause and is NOT a supplement table, so the variants live here in data exploration.
# Each variant is c_07 with one line changed, made by text substitution on the script itself so
# there is exactly one copy of the model code. Variants:
#
#   bird window   21-28, 14-28, 8-28 d (the reported model uses 15-28 d)
#   weather 0-7   weather term restricted to the first week, birds 7-28 d
#   drop-one      runoff dropped; precipitation dropped - the two water terms carry opposite signs,
#                 and each is refitted without the other to show that is not shared variance
#
# Each writes si_interaction_primary_<tag>.RDS alongside c_07's si_interaction_primary.RDS.

project.folder = paste0(print(here::here()), '/')
src <- readLines(paste0(project.folder, 'create_folder_structure.R'))
eval(parse(text = paste(src[!grepl("^\\s*rm\\(list", src)], collapse = "\n")))
BASE <- readLines(paste0(project.folder, "02_code/2c_models/c_07_interaction.R"))

sub_line <- function(txt, from, to) {
  hit <- grepl(from, txt, fixed = TRUE)
  if (!any(hit)) stop("no line matches: ", from)
  sub(from, to, txt, fixed = TRUE)
}
VARIANTS <- list(
  `w21-28`           = list(c("BIRDLAG <- 15:28", "BIRDLAG <- 21:28")),
  `w14-28`           = list(c("BIRDLAG <- 15:28", "BIRDLAG <- 14:28")),
  `w8-28`            = list(c("BIRDLAG <- 15:28", "BIRDLAG <- 8:28")),
  `met0-7_bird7-28`  = list(c("BIRDLAG <- 15:28", "BIRDLAG <- 7:28"),
                            c("MAXLAG  <- 14", "MAXLAG  <- 7"),
                            c("WINDOW  <- list(`0-7 days` = 0:7, `8-14 days` = 8:14, `0-14 days` = 0:14)",
                              "WINDOW  <- list(`0-7 days` = 0:7)"),
                            c('ARGLAG  <- list(fun = "strata", breaks = 8)', 'ARGLAG  <- list(fun = "strata", df = 1)')),
  `norunoff`         = list(c("VARS    <- INLA_PRIMARY_MET", 'VARS    <- setdiff(INLA_PRIMARY_MET, "runoff")')),
  `noprecip`         = list(c("VARS    <- INLA_PRIMARY_MET", 'VARS    <- setdiff(INLA_PRIMARY_MET, "precipitation")')))

for (tag in names(VARIANTS)) {
  txt <- BASE
  for (s in VARIANTS[[tag]]) txt <- sub_line(txt, s[1], s[2])
  txt <- sub_line(txt, 'saveRDS(res, paste0(objects_folder, "si_interaction_primary.RDS"))',
                  sprintf('saveRDS(res, paste0(objects_folder, "si_interaction_primary_%s.RDS"))', tag))
  # the base script starts by printing here() and sourcing; running it in a fresh session keeps
  # each variant's globals from leaking into the next
  f <- tempfile(fileext = ".R"); writeLines(txt, f)
  cat(sprintf("\n=== variant %s ===\n", tag)); flush.console()
  status <- system2("Rscript", f, stdout = "", stderr = "")
  if (status != 0) warning("variant ", tag, " failed")
}
cat("\nwrote si_interaction_primary_<tag>.RDS for:", paste(names(VARIANTS), collapse = ", "), "\n")
