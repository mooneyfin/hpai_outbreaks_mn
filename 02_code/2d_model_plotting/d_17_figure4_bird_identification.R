############ FIGURE 4 — THE LAG DISTRIBUTION IS NOT IDENTIFIED ############
# Both bird terms, both lag windows, cumulative within disjoint weekly windows.
#
# The point of the figure is the contradiction. Same 3,870 rows, same 167 strata, same design,
# prior and met adjustment - only the bird lag window changes - and the waterfowl gradient
# inverts. Rising 0.78 -> 2.53 at 0-28 d, falling 1.83 -> 1.42 at 0-56 d, with credible
# intervals at lag 0-7 that don't overlap. Showing one window and not the other would be
# choosing the story; showing both says what the data can actually support, which is THAT
# waterfowl abundance is associated with spillover, not WHEN.
#
# Predators is here as the comparison, not as a covariate. It was dropped from the primary
# because scavengers feeding on infected carcasses are a plausible mediator, so adjusting for it
# would be over-adjustment. Each bird term therefore gets its OWN model alongside the met terms
# - two total effects, compared like for like. Within-stratum correlation between them is 0.34,
# so this isn't a collinearity dodge, it's the mediator argument.
#
# Cumulative, not lag-specific: it's the convention used everywhere else in the paper, and the
# per-day shape is the unstable thing this figure exists to document.
#
# Caption needs the attenuation note: case and referent share more exposure history the further
# out you go (55% at 28 d, 77% at 56 d), so late windows are pulled toward the null. That makes
# a rising late estimate credible and a flat one ambiguous - it does not rescue a timing claim.

rm(list = ls())
project.folder = paste0(print(here::here()), '/')
source(paste0(project.folder, 'create_folder_structure.R'))
source(paste0(functions.folder, 'script_initiate.R'))
suppressMessages({library(INLA); library(Matrix); library(viridisLite)})
source(paste0(functions.folder, 'inla_dlnm_helpers.R'))

# same prior and exposure set as the primary model; both come from inla_dlnm_helpers.R
PREC <- INLA_PREC_PRIMARY; NS2 <- list(fun = "ns", df = 2); LIN <- list(fun = "lin")
SH   <- INLA_PRIMARY_SHOCK
BIRDS <- c(anseriformes = "Anseriformes", predators = "Predators")
WIN  <- list(c(0,7), c(8,14), c(15,21), c(22,28), c(29,35), c(36,42), c(43,49), c(50,56))
WLAB <- vapply(WIN, function(w) sprintf("%d-%d", w[1], w[2]), character(1))
MODLAB <- c(`28` = "Bird lag 0-28 days", `56` = "Bird lag 0-56 days")

PANEL_ORD <- INLA_PANEL_ORD
OK  <- setNames(viridisLite::viridis(3, begin = 0.05, end = 0.78), PANEL_ORD)
SHP <- setNames(c(16, 17, 15), PANEL_ORD)

D <- setDT(readRDS(paste0(objects_folder, "case_crossover_df_shock_post7.RDS")))
if (anyDuplicated(names(D))) D <- D[, which(!duplicated(names(D))), with = FALSE]
sub <- function(pn) {
  if (pn == "Minnesota") D[state == "Minnesota"]
  else if (pn == INLA_PANEL_LAB[["Pooled"]]) copy(D)
  else D[state != "Minnesota"]
}

RDS <- paste0(objects_folder, "figure3_bird_curves.RDS")   # v2: six exposures, primary prior, lag-specific added
REFIT <- !file.exists(RDS)

grid <- CJ(bird = names(BIRDS), ml = c(28L, 56L), panel = PANEL_ORD, sorted = FALSE)
res <- if (!REFIT) setDT(readRDS(RDS)) else rbindlist(lapply(seq_len(nrow(grid)), function(i) {
  bv <- grid$bird[i]; ml <- grid$ml[i]; pn <- grid$panel[i]
  d  <- sub(pn)
  av <- c(setNames(rep(list(NS2), length(SH)), SH), setNames(list(LIN), bv))
  o  <- assemble_cc_dlnm(d, met_vars = c(SH, bv), weekly_vars = character(0),
                         max_lag = ml, max_lag_by_var = setNames(as.list(rep(28, length(SH))), SH),
                         argvar_by_var = av, arglag = list(fun = "ns", df = 2))
  f  <- suppressWarnings(fit_cc_inla_dlnm(o, fixed_prec = PREC))
  cat(sprintf("  %-12s lag 0-%2d  %-20s n=%3d  WAIC %.1f\n", BIRDS[[bv]], ml, pn,
              sum(d$outbreak_binary), f$waic$waic))
  # two curves per fit: cumulative from lag 0 to L (start-point dependent by construction) and
  # the LAG-SPECIFIC per-day contribution at L, which has no start point and shows where in the
  # lag the effect sits
  rbindlist(lapply(0:ml, function(L) {
    e <- met_lag_effect(f, o, bv, 0:L); s <- met_lag_effect(f, o, bv, L)
    data.table(bird = BIRDS[[bv]], model = MODLAB[[as.character(ml)]], panel = pn,
               lag = L, waic = f$waic$waic,
               RR = e[["OR"]], low = e[["low"]], high = e[["high"]],
               RRl = s[["OR"]], lowl = s[["low"]], highl = s[["high"]])
  }))
}))
res[, panel := norm_panel(panel)]
res[, `:=`(panel = factor(panel, levels = PANEL_ORD),
           bird  = factor(bird,  levels = unname(BIRDS)),
           model = factor(model, levels = unname(MODLAB)))]
stopifnot(!anyNA(res$panel))
if (REFIT) saveRDS(res, RDS)

# ---- 3. Figures ----
# LAG-SPECIFIC curves, not cumulative-from-zero. A cumulative curve's starting point depends on
# where the fitted lag range starts, so the 0-28 and 0-56 d fits looked as if they disagreed about
# WHEN the effect begins; that was an artefact of the summary. The per-day rate ratio has no start
# point: it shows where in the lag the waterfowl effect sits, and whether it keeps rising past 28 d.
# The cumulative curves stay in the RDS (RR/low/high) for anyone who wants them.
lab_num <- scales::label_number(drop0trailing = TRUE)
lagplot <- function(d, colour_by, title = NULL) {
  ggplot(d, aes(lag, RRl, colour = .data[[colour_by]], fill = .data[[colour_by]])) +
    geom_hline(yintercept = 1, colour = "grey45", linewidth = 0.3) +
    geom_ribbon(aes(ymin = lowl, ymax = highl), alpha = 0.14, colour = NA) +
    geom_line(linewidth = 0.85) +
    scale_x_continuous(breaks = seq(0, 56, 7)) +
    scale_y_continuous(trans = "log", breaks = c(0.8, 0.9, 1, 1.1, 1.2, 1.5), labels = lab_num) +
    labs(x = "Lag (days)", y = "Lag-specific rate ratio per +0.5 SD bird abundance",
         colour = title, fill = title) +
    theme_classic(base_family = "Avenir") +
    theme(legend.position = "top", legend.justification = "left",
          legend.title = element_text(family = "Avenir Heavy", size = 13, vjust = 0.5),
          legend.text = element_text(size = 11),
          strip.text = element_text(family = "Avenir Heavy", size = 11), strip.background = element_blank(),
          axis.title = element_text(size = 11), axis.text = element_text(size = 9.5, colour = "black"),
          panel.border = element_rect(fill = NA, colour = "grey20", linewidth = 0.4))
}
# Figure 3 (main text): Minnesota, the 0-28 d fit that matches the primary model, both bird terms
mn <- res[panel == "Minnesota" & model == MODLAB[["28"]]]
f3 <- lagplot(mn, "bird", "Minnesota, 90-day shock model, lag 0-28 d") +
  scale_colour_manual(values = c(Anseriformes = "#440154", Predators = "#5DC863")) +
  scale_fill_manual(values = c(Anseriformes = "#440154", Predators = "#5DC863"))
ggsave_spark(file.path(figures_main_folder, "figure3_bird_lag_response.png"), f3, width = 8, height = 5)
# Figure S4: all panels, both lag ranges - does the effect keep rising past 28 d anywhere?
fS4 <- lagplot(res, "panel", NULL) + facet_grid(model ~ bird, scales = "free_x") +
  scale_colour_manual(values = OK) + scale_fill_manual(values = OK)
ggsave_spark(file.path(figures_main_folder, "figureS4_bird_lag_response_all_panels.png"), fS4, width = 12, height = 8)
cat("wrote figure3_bird_lag_response.png (Minnesota) and figureS4_bird_lag_response_all_panels.png\n")
res[, cell := sprintf("%.2f (%.2f, %.2f)", RRl, lowl, highl)]
print(dcast(res[bird == "Anseriformes" & lag %in% c(0, 7, 14, 21, 28, 42, 56)], lag ~ model + panel, value.var = "cell"))
