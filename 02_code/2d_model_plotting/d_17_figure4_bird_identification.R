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

PREC <- 4; NS2 <- list(fun = "ns", df = 2); LIN <- list(fun = "lin")
SH   <- c("temperature_shock", "precipitation_shock", "soil_moisture_shock", "wind_speed_shock")
BIRDS <- c(anseriformes = "Anseriformes", predators = "Predators")
WIN  <- list(c(0,7), c(8,14), c(15,21), c(22,28), c(29,35), c(36,42), c(43,49), c(50,56))
WLAB <- vapply(WIN, function(w) sprintf("%d-%d", w[1], w[2]), character(1))
MODLAB <- c(`28` = "Bird lag 0-28 days", `56` = "Bird lag 0-56 days")

PANEL_ORD <- c("Minnesota", "Other northern Mississippi Flyway states", "Pooled")
OK  <- setNames(viridisLite::viridis(3, begin = 0.05, end = 0.78), PANEL_ORD)
SHP <- setNames(c(16, 17, 15), PANEL_ORD)

D <- setDT(readRDS(paste0(objects_folder, "case_crossover_df_shock_post7.RDS")))
if (anyDuplicated(names(D))) D <- D[, which(!duplicated(names(D))), with = FALSE]
sub <- function(pn) {
  if (pn == "Minnesota") D[state == "Minnesota"]
  else if (pn == "Pooled") copy(D)
  else D[state != "Minnesota"]
}

RDS <- paste0(objects_folder, "figure4_bird_identification.RDS")
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
  # cumulative from lag 0 up to L, for every L - this is the curve, x is days
  rbindlist(lapply(0:ml, function(L) {
    e <- met_lag_effect(f, o, bv, 0:L)
    data.table(bird = BIRDS[[bv]], model = MODLAB[[as.character(ml)]], panel = pn,
               lag = L, waic = f$waic$waic,
               RR = e[["OR"]], low = e[["low"]], high = e[["high"]])
  }))
}))
res[panel == "Other flyway states", panel := "Other northern Mississippi Flyway states"]
res[, `:=`(panel = factor(panel, levels = PANEL_ORD),
           bird  = factor(bird,  levels = unname(BIRDS)),
           model = factor(model, levels = unname(MODLAB)))]
stopifnot(!anyNA(res$panel))
if (REFIT) saveRDS(res, RDS)

# BOTH lag windows are plotted, which is the whole point of the figure per the header above.
# The earlier version showed only the 0-56 d fit on the grounds that the cumulative curves agree
# at lag 28 - and pooled they do, 3.16 (1.95, 5.12) against 3.41 (2.10, 5.55). But the figure's
# visual claim is the SHAPE, not the endpoint, and the shapes do not agree in Minnesota:
#
#   lag 0    0-28 d model  0.94 (0.87, 1.02)   0-56 d model  1.08 (1.04, 1.12)
#   lag 7                  0.78 (0.51, 1.20)                 1.83 (1.39, 2.42)
#
# Neither pair overlaps. The 0-56 d curve excludes 1 at every lag; the 0-28 d curve - the window
# the PRIMARY results are read from - crosses 1 until lag 24. Showing only the first would let a
# reader take a smooth monotone rise from lag 0 as established when it is an artefact of the lag
# range. Pooled is still dropped: it's a weighted blend of the other two, not a third thing.
plotdat <- res[panel != "Pooled"]
plotdat[, panel := droplevels(panel)]

lab_num <- scales::label_number(drop0trailing = TRUE)
fig4 <- ggplot(plotdat, aes(lag, RR, colour = panel, fill = panel)) +
  geom_hline(yintercept = 1, colour = "grey45", linewidth = 0.3) +
  geom_ribbon(aes(ymin = low, ymax = high), alpha = 0.12, colour = NA) +
  geom_line(aes(linetype = panel), linewidth = 0.85) +
  facet_grid(model ~ bird) +
  scale_colour_manual(values = OK) + scale_fill_manual(values = OK) +
  scale_linetype_manual(values = setNames(c("solid", "dashed", "dotdash"), PANEL_ORD),
                        drop = TRUE) +
  scale_x_continuous(breaks = seq(0, 56, 7)) +
  # log y: the cumulative climbs several-fold, and a linear axis would flatten the early lags
  # into nothing while stretching the tail
  scale_y_continuous(trans = "log", breaks = c(0.5, 1, 2, 5, 10, 20), labels = lab_num) +
  labs(x = "Lag (days)",
       y = "Cumulative rate ratio per +0.5 SD bird abundance",
       colour = NULL, fill = NULL, linetype = NULL) +
  theme_classic(base_family = "Avenir") +
  # legend tucked into the top right of the Predators panel, which is empty - that curve sits
  # flat on the null the whole way
  # legend moves to the top: with a 2x2 grid there is no longer an empty panel to tuck it into
  theme(legend.position = "top", legend.justification = "left",
        strip.text = element_text(family = "Avenir Heavy", size = 11),
        axis.title = element_text(size = 12),
        axis.text = element_text(size = 9.5, colour = "black"),
        legend.text = element_text(size = 14),
        legend.key.size = unit(1.6, "lines"),
        panel.border = element_rect(fill = NA, colour = "grey20", linewidth = 0.4))
ggsave_spark(file.path(figures_main_folder, "figure4_bird_cumulative_curves.png"),
             fig4, width = 12, height = 9)

res[, cell := sprintf("%.2f (%.2f, %.2f)", RR, low, high)]
for (bv in unname(BIRDS)) {
  cat(sprintf("\n=== %s : cumulative 0 to L ===\n", bv))
  print(dcast(res[bird == bv & lag %in% c(7, 14, 21, 28, 42, 56)],
              lag ~ model + panel, value.var = "cell"))
}
cat("\nwrote figure4_bird_cumulative_curves.png\n")
