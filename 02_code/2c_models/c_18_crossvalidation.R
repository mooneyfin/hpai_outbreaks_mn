############ SUPPLEMENTARY TABLE S8 — CROSS-VALIDATION ############
# Stratum-level 5-fold CV under the final specification. The score is the conditional
# log-likelihood of the observed case day within each held-out stratum, which is what a
# case-crossover actually predicts; the stratum intercept cancels, so it never has to be
# estimated for a stratum the model has not seen.
#
# Comparisons are PAIRED on the stratum, because the fold assignment is shared. Without the
# pairing there is no way to tell a real difference from fold noise at this sample size.
#
# Note the one comparison that is NOT valid and is deliberately absent: designs with
# different washouts cannot be compared this way. The conditional log-likelihood scales with
# the number of referents per stratum, so a design with fewer referents scores better for
# arithmetic reasons rather than predictive ones.

rm(list = ls())
project.folder = paste0(print(here::here()), '/')
source(paste0(project.folder, 'create_folder_structure.R'))
source(paste0(functions.folder, 'script_initiate.R'))
suppressMessages({library(INLA); library(Matrix); library(flextable); library(officer)})
INLA::inla.setOption(num.threads = 1)
source(paste0(functions.folder, 'inla_dlnm_helpers.R'))

SEED <- 20240517; K <- 5; PREC <- 4
NS2  <- list(fun = "ns", df = 2); LIN <- list(fun = "lin")
SH   <- c("temperature_shock", "precipitation_shock", "soil_moisture_shock", "wind_speed_shock")
LV   <- c("temperature", "precipitation", "soil_moisture", "wind_speed")

D <- setDT(readRDS(paste0(objects_folder, "case_crossover_df_shock_post7.RDS")))
if (anyDuplicated(names(D))) D <- D[, which(!duplicated(names(D))), with = FALSE]

# linear predictor from the exposure terms only; the stratum intercept cancels in the
# conditional likelihood so it is deliberately excluded
eta_fixed <- function(fit, obj) {
  cbn <- unlist(obj$cb_cols)
  b <- fit$summary.fixed[match(cbn, rownames(fit$summary.fixed)), "mean"]
  as.numeric(as.matrix(obj$data[, cbn, drop = FALSE]) %*% b)
}
cond_ll <- function(eta, y, st) {
  d <- data.table(eta, y, st)
  d[, .(ll = eta[y == 1][1] - (max(eta) + log(sum(exp(eta - max(eta)))))), by = st]
}

SPECS <- list(
  "Anomalies + waterfowl (primary)" = c(SH, "anseriformes"),
  "Absolute conditions + waterfowl" = c(LV, "anseriformes"),
  "Anomalies only (no waterfowl)"   = SH,
  "Waterfowl only"                  = "anseriformes",
  "Soil moisture anomaly + waterfowl" = c("soil_moisture_shock", "anseriformes"))

run_cv <- function(dat, label) {
  set.seed(SEED)
  st   <- unique(dat$stratum_id)
  fold <- setNames(sample(rep_len(1:K, length(st))), st)
  rbindlist(lapply(names(SPECS), function(nm) {
    mv <- SPECS[[nm]]
    av <- setNames(lapply(mv, function(v) if (v == "anseriformes") LIN else NS2), mv)
    ll <- rbindlist(lapply(1:K, function(k) {
      tr <- dat[stratum_id %in% names(fold)[fold != k]]
      te <- dat[stratum_id %in% names(fold)[fold == k]]
      o_tr <- assemble_cc_dlnm(tr, met_vars = mv, weekly_vars = character(0), max_lag = 28,
                               argvar_by_var = av, arglag = list(fun = "ns", df = 2))
      o_te <- assemble_cc_dlnm(te, met_vars = mv, weekly_vars = character(0), max_lag = 28,
                               argvar_by_var = av, arglag = list(fun = "ns", df = 2))
      f_tr <- suppressWarnings(fit_cc_inla_dlnm(o_tr, fixed_prec = PREC))
      cond_ll(eta_fixed(f_tr, o_te), o_te$data$outbreak_binary, te$stratum_id)
    }))
    cat(sprintf("  %-36s %s  CV logLik %.1f\n", nm, label, sum(ll$ll)))
    data.table(panel = label, spec = nm, stratum = ll$st, ll = ll$ll)
  }))
}

cat("running 5-fold stratum CV\n")
cv <- rbindlist(lapply(c("Minnesota", "Pooled"), function(p)
  run_cv(if (p == "Minnesota") D[state == "Minnesota"] else copy(D), p)))

# an equal-odds model, for scale: it predicts every day in a stratum equally, so its
# per-stratum log-likelihood is -log(number of days in that stratum). Count those from the
# panel itself - counting rows of `cv` gives 1 per stratum and a meaningless zero.
nday <- rbind(
  D[state == "Minnesota", .(nd = .N), by = stratum_id][, panel := "Minnesota"],
  D[, .(nd = .N), by = stratum_id][, panel := "Pooled"])
base <- nday[, .(nul = sum(-log(nd))), by = panel]

REF <- names(SPECS)[1]
tab <- rbindlist(lapply(c("Minnesota", "Pooled"), function(p) {
  x <- cv[panel == p]
  r <- x[spec == REF, .(stratum, ref = ll)]
  m <- merge(x, r, by = "stratum")
  s <- m[, .(cv = sum(ll), d = mean(ll - ref),
             z = mean(ll - ref) / (sd(ll - ref) / sqrt(.N))), by = spec]
  s[, .(Predictor = paste0(spec, " - ", p),
        `CV log-likelihood` = sprintf("%.1f", cv),
        `vs equal-odds` = sprintf("%+.1f", cv - base[panel == p]$nul),
        `Paired difference vs primary` = ifelse(spec == REF, "reference", sprintf("%+.4f", d)),
        `z` = ifelse(spec == REF, "-", sprintf("%+.2f", z)))]
}))
saveRDS(list(summary = tab, per_stratum = cv),
        paste0(objects_folder, "tableS8_crossvalidation.RDS"))

ft <- flextable(tab) |>
  add_header_row(values = c("Specification", "5-fold stratum cross-validation"),
                 colwidths = c(1, ncol(tab) - 1)) |>
  fontsize(size = 11, part = "all") |> autofit() |>
  align(align = "center", part = "all") |> theme_booktabs()
read_docx() |> body_add_flextable(ft) |>
  print(target = file.path(tables_main_folder, "TableS8.docx"))

cat("\n=== Table S5 ===\n"); print(tab)
cat("\nwrote TableS8.docx\n")
