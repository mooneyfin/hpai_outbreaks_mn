############ TABLE 2 AND SUPPLEMENTARY TABLES S2-S4 ############
# Table 2  - primary results, cumulative ORs, all three panels
# Table S2 - prior sensitivity + prior-predictive check (justifies N(0, 0.5^2))
# Table S3 - level-parameterisation results, the companion to Figure 3a
# Table S4 - referent design comparison, time-stratified designs only
#
# Bare tables: no titles, no footnotes. Captions live in the manuscript text.

rm(list = ls())
project.folder = paste0(print(here::here()), '/')
source(paste0(project.folder, 'create_folder_structure.R'))
source(paste0(functions.folder, 'script_initiate.R'))
suppressMessages({library(INLA); library(Matrix); library(flextable); library(officer)})
INLA::inla.setOption(num.threads = 1)
source(paste0(functions.folder, 'inla_dlnm_helpers.R'))

fits   <- readRDS(paste0(objects_folder, "primary_shock_level_fits.RDS"))
PANELS <- c("Minnesota", "Other flyway states", "Pooled")
KEYS   <- c(Minnesota = "Minnesota", `Other flyway states` = "Other flyway",
            Pooled = "Pooled")
WCOLS  <- names(INLA_LAG_WINDOWS)

# bold anything whose interval excludes 1
is_sig <- function(x) {
  m <- regmatches(x, regexec("\\(([0-9.]+), ([0-9.]+)\\)", x))
  vapply(m, function(z) length(z) == 3 && (as.numeric(z[2]) > 1 | as.numeric(z[3]) < 1), TRUE)
}
ms_ft <- function(t, span) {
  ft <- flextable(t) |>
    add_header_row(values = c("Predictor", span), colwidths = c(1, ncol(t) - 1)) |>
    fontsize(size = 11, part = "all") |> autofit() |>
    align(align = "center", part = "all") |> theme_booktabs()
  for (j in setdiff(names(t), "Predictor")) {
    hit <- which(is_sig(t[[j]]))
    if (length(hit)) ft <- bold(ft, i = hit, j = j)
  }
  ft
}
wr <- function(ft, f) read_docx() |> body_add_flextable(ft) |>
  print(target = file.path(tables_main_folder, f))

# 1. Table 2 - primary model, cumulative AND lag-specific. The two are deterministically
#    linked (a window's cumulative is the sum of its lag-specific terms), but the convention
#    in a DLNM paper is to show both: the per-day contribution and what it accumulates to.
ACOLS <- paste0("Lag ", INLA_LAG_ANCHORS, " days")
t2 <- rbindlist(lapply(PANELS, function(p) {
  cum <- copy(fits[[paste0(KEYS[[p]], "|shock")]]$cumulative)
  lag <- copy(fits[[paste0(KEYS[[p]], "|shock")]]$lagspecific)
  rbind(cbind(data.table(Panel = p, Estimate = "Cumulative"),
              setnames(cum, WCOLS, paste0("W", seq_along(WCOLS)))),
        cbind(data.table(Panel = p, Estimate = "Lag-specific"),
              setnames(lag[, c("Predictor", ACOLS[1:4]), with = FALSE],
                       ACOLS[1:4], paste0("W", 1:4))))
}))
setnames(t2, paste0("W", 1:4), WCOLS)
setcolorder(t2, c("Panel", "Estimate", "Predictor", WCOLS))
saveRDS(t2, paste0(objects_folder, "table2_primary.RDS"))
wr(ms_ft(t2, "Odds ratio (95% credible interval) per +0.5 SD"),
   "Table2_primary_results.docx")

# 2. Table S3 - level model, same layout
s3 <- rbindlist(lapply(PANELS, function(p) {
  d <- copy(fits[[paste0(KEYS[[p]], "|level")]]$cumulative)
  cbind(data.table(Panel = p), d)
}))
setcolorder(s3, c("Panel", "Predictor", WCOLS))
saveRDS(s3, paste0(objects_folder, "tableS3_level_model.RDS"))
wr(ms_ft(s3, "Cumulative odds ratio (95% credible interval) per +0.5 SD"),
   "TableS3_level_model.docx")

# 3. Table S2 - prior sensitivity. WAIC has an interior optimum at prec 64, but that optimum
#    shrinks the effects; the prior is chosen on the prior-predictive column instead.
D  <- setDT(readRDS(paste0(objects_folder, "case_crossover_df_shock_post7.RDS")))
if (anyDuplicated(names(D))) D <- D[, which(!duplicated(names(D))), with = FALSE]
SH <- c("temperature_shock", "precipitation_shock", "soil_moisture_shock",
        "wind_speed_shock")
AV <- c(setNames(rep(list(list(fun = "ns", df = 2)), length(SH)), SH),
        list(anseriformes = list(fun = "lin")))
o  <- assemble_cc_dlnm(D, met_vars = c(SH, "anseriformes"), weekly_vars = character(0),
                       max_lag = 28, argvar_by_var = AV, arglag = list(fun = "ns", df = 2))
GRID <- c(0.001, 1, 4, 16, 64, 256)
set.seed(1)
cb <- o$bases[["soil_moisture_shock"]]
ob <- function(x, a) do.call(dlnm::onebasis, c(list(x = x), a))
v  <- as.numeric(ob(0.5, attr(cb, "argvar")) - ob(0, attr(cb, "argvar")))
W  <- colSums(ob(0:7, attr(cb, "arglag")), na.rm = TRUE)
cvec <- as.vector(t(outer(v, W)))

s2 <- rbindlist(lapply(GRID, function(pr) {
  f <- suppressWarnings(fit_cc_inla_dlnm(o, fixed_prec = pr))
  e <- met_lag_effect(f, o, "soil_moisture_shock", 0:7)
  b <- met_lag_effect(f, o, "anseriformes", 22:28)
  B  <- matrix(rnorm(20000 * length(cvec), 0, 1 / sqrt(pr)), ncol = length(cvec))
  OR <- exp(as.numeric(B %*% cvec))
  data.table(Predictor = sprintf("N(0, %.3g^2)", 1 / sqrt(pr)),
             WAIC = sprintf("%.1f", f$waic$waic),
             `Soil moisture 0-7 days` = sprintf("%.2f (%.2f, %.2f)", e[["OR"]], e[["low"]], e[["high"]]),
             `Anseriformes 22-28 days` = sprintf("%.2f (%.2f, %.2f)", b[["OR"]], b[["low"]], b[["high"]]),
             `Prior-implied OR, 95% range` = sprintf("%.2f to %.2f",
                                                     quantile(OR, .025), quantile(OR, .975)))
}))
saveRDS(s2, paste0(objects_folder, "tableS2_prior_sensitivity.RDS"))
wr(ms_ft(s2, "Prior on the crossbasis coefficients"), "TableS2_prior_sensitivity.docx")

# 4. Table S4 - referent design, time-stratified only
DES <- c("No washout" = "case_crossover_df_timestrat_month.RDS",
         "Post-only 7-day post exclusion" = "case_crossover_df_timestrat_month_post7.RDS",
         "Post-only 14-day post exclusion" = "case_crossover_df_timestrat_month_post14.RDS")
ex <- setDT(readRDS(paste0(objects_folder, "daily_shock_predictors.RDS")))
ex[, date := as.Date(date)]
SHOCK <- setdiff(names(ex), c("zone_id", "date"))
attach_shocks <- function(f) {
  cc <- setDT(readRDS(paste0(objects_folder, f)))
  if (anyDuplicated(names(cc))) cc <- cc[, which(!duplicated(names(cc))), with = FALSE]
  cc[, date := as.Date(date)]
  for (L in 0:28) {
    src <- ex[, c(.(zone_id = zone_id, date = date + L), .SD), .SDcols = SHOCK]
    setnames(src, SHOCK, paste0(SHOCK, "_Lag", L))
    cc <- merge(cc, src, by = c("zone_id", "date"), all.x = TRUE, sort = FALSE)
  }
  for (vv in SHOCK) for (L in 0:28) {
    cn <- paste0(vv, "_Lag", L); set(cc, j = cn, value = as.numeric(scale(cc[[cn]])))
  }
  need <- paste0(rep(SHOCK, each = 29), "_Lag", 0:28)
  cc <- cc[complete.cases(cc[, ..need])]
  cc[, `:=`(nc = sum(outbreak_binary), nr = sum(outbreak_binary == 0)), by = stratum_id]
  cc[nc == 1L & nr > 0L]
}
g <- function(t, r, c) { x <- t[Predictor == r][[c]]; if (length(x)) x else NA_character_ }
s4 <- rbindlist(lapply(names(DES), function(dn) {
  P <- attach_shocks(DES[[dn]])
  ct  <- as.Date(sub(".*_", "", P$stratum_id), "%Y%m%d")
  bal <- mean(as.integer(as.Date(P$date) - ct)[P$outbreak_binary == 0])
  rbindlist(lapply(c("Minnesota", "Pooled"), function(sn) {
    d <- if (sn == "Minnesota") P[state == "Minnesota"] else copy(P)
    oo <- assemble_cc_dlnm(d, met_vars = c(SH, "anseriformes"), weekly_vars = character(0),
                           max_lag = 28, argvar_by_var = AV, arglag = list(fun = "ns", df = 2))
    ff <- suppressWarnings(fit_cc_inla_dlnm(oo, fixed_prec = 4))
    tt <- lag_window_table(ff, oo, "cumulative", labels = fits[["Pooled|shock"]]$labels)
    data.table(Predictor = paste0(dn, " - ", sn),
               `Referents per case` = sprintf("%.1f", sum(d$outbreak_binary == 0) / sum(d$outbreak_binary)),
               `Referent balance (days)` = sprintf("%+.1f", bal),
               `Soil moisture 0-7 days` = g(tt, "Soil moisture", WCOLS[1]),
               `Anseriformes 22-28 days` = g(tt, "Anseriformes", WCOLS[4]))
  }))
}))
saveRDS(s4, paste0(objects_folder, "tableS4_referent_design.RDS"))
wr(ms_ft(s4, "Time-stratified referent designs"), "TableS4_referent_design.docx")

cat("\n=== Table 2 ===\n");  print(t2)
cat("\n=== Table S2 ===\n"); print(s2)
cat("\n=== Table S4 ===\n"); print(s4)
cat("\nwrote Table2, TableS2, TableS3, TableS4 (.docx + .RDS)\n")
