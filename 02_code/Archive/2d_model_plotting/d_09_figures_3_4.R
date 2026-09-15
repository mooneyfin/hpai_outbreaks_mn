############ FIGURE 3 ############
# Fig 3: cumulative exposure-response, (a) absolute conditions and (b) anomalies. The point of
#        putting them side by side is that it answers "why shocks?" in one look.
# Fig 4: the waterfowl lag profile across the three panels. This is the headline result.
# Both read the cached fits from c_17.

rm(list = ls())
project.folder = paste0(print(here::here()), '/')
source(paste0(project.folder, 'create_folder_structure.R'))
source(paste0(functions.folder, 'script_initiate.R'))
suppressMessages({library(INLA); library(Matrix); library(patchwork); library(viridisLite)})
source(paste0(functions.folder, 'inla_dlnm_helpers.R'))

fits <- readRDS(paste0(objects_folder, "primary_shock_level_fits.RDS"))
lab_num <- scales::label_number(drop0trailing = TRUE)   # 1 and 1.5, never 1.0
# viridis, dark end first so Minnesota reads as the focal panel
PANEL_ORD <- c("Minnesota", "Other flyway states", "Pooled")
OK <- setNames(viridisLite::viridis(3, begin = 0.05, end = 0.78), PANEL_ORD)
SHP <- setNames(c(16, 17, 15), PANEL_ORD)
RELABEL <- c(Minnesota = "Minnesota", `Other flyway` = "Other flyway states",
             Pooled = "Pooled")

# 1. pull cumulative estimates out of a fit as tidy numbers rather than the formatted strings
cum_tidy <- function(r, panel, kind) {
  rbindlist(lapply(names(INLA_LAG_WINDOWS), function(k) {
    w <- INLA_LAG_WINDOWS[[k]]
    rbindlist(lapply(c(r$obj$met_vars), function(v) {
      e <- met_lag_effect(r$fit, r$obj, v, seq(w[1], w[2]))
      data.table(panel = panel, kind = kind, variable = unname(r$labels[v]),
                 window = k, OR = e[["OR"]], low = e[["low"]], high = e[["high"]])
    }))
  }))
}
lagspec_tidy <- function(r, panel, v = "anseriformes") {
  rbindlist(lapply(0:r$obj$max_lag, function(a) {
    e <- met_lag_effect(r$fit, r$obj, v, a)
    data.table(panel = panel, lag = a, OR = e[["OR"]], low = e[["low"]], high = e[["high"]])
  }))
}

# 2. Figure 3 - levels vs shocks, all three panels. This absorbs what used to be a separate
#    flyway forest plot: that figure's pooled column was just this figure's row (b), so the
#    two were showing the same 20 estimates twice.
both <- rbindlist(lapply(names(RELABEL), function(p)
  rbind(cum_tidy(fits[[paste0(p, "|level")]], RELABEL[[p]], "level"),
        cum_tidy(fits[[paste0(p, "|shock")]], RELABEL[[p]], "shock"))))
# strip the " shock" suffix so both rows share y labels and read as a like-for-like comparison
both[, variable := sub(" shock$", "", variable)]
both[, variable := sub("^Surface runoff$", "Runoff", variable)]
ORD <- c("Temperature", "Precipitation", "Soil moisture", "Wind speed", "Anseriformes")
both <- both[variable %in% ORD]
both[, variable := factor(variable, levels = rev(ORD))]
both[, window := factor(window, levels = names(INLA_LAG_WINDOWS))]
both[, panel := factor(panel, levels = PANEL_ORD)]

# wide enough to hold every interval (lowest bound 0.45, highest 3.27) so nothing is clipped
XLIM <- c(0.5, 3.5)
both[, `:=`(low_draw = low, high_draw = high, clipped_hi = FALSE)]

panel_plot <- function(d, ttl) {
  ggplot(d, aes(OR, variable, colour = panel, shape = panel)) +
    geom_vline(xintercept = 1, colour = "grey45", linewidth = 0.3) +
    # reverse = TRUE so Minnesota sits at the top of each group and Pooled at the bottom,
    # matching the legend order rather than flipping it
    geom_errorbarh(aes(xmin = low_draw, xmax = high_draw), height = 0, linewidth = 0.55,
                   position = position_dodge(width = 0.66, reverse = TRUE)) +
    geom_point(size = 2.1, position = position_dodge(width = 0.66, reverse = TRUE)) +
    facet_wrap(~ window, nrow = 1) +
    scale_colour_manual(values = OK) +
    scale_shape_manual(values = SHP) +
    # log x: ratios are multiplicative, so 0.5 and 2.0 are the same effect in opposite
    # directions. on a linear axis they sit at different distances from the null and the
    # harmful side looks bigger than it is.
    scale_x_continuous(trans = "log", breaks = c(0.5, 1, 2, 3), labels = lab_num) +
    # coord_cartesian, NOT scale limits: `limits` DROPS any estimate whose interval leaves the
    # range, which silently deleted four whiskers here (MN temperature down to 0.46, MN
    # Anseriformes up to 3.66 and 4.00). This clips the drawing instead, so a bar that runs
    # past the edge reads as running past the edge.
    coord_cartesian(xlim = XLIM) +
    labs(x = "Rate ratio per +0.5 standard deviation increase in exposure",
         y = NULL, colour = ttl, shape = ttl) +
    theme_classic(base_family = "Avenir") +
    theme(legend.title = element_text(family = "Avenir Heavy", size = 14,
                                      vjust = 0.5, margin = margin(r = 14)),
          strip.text = element_text(family = "Avenir Heavy", size = 11),
          strip.background = element_blank(),
          axis.text = element_text(size = 10, colour = "black"),
          axis.title = element_text(size = 12),
          legend.position = "top", legend.justification = "left",
          legend.text = element_text(size = 14),
          legend.key.size = unit(1.4, "lines"),
          legend.box.spacing = unit(2, "pt"),
          panel.border = element_rect(fill = NA, colour = "grey20", linewidth = 0.4))
}
# both rows carry their own legend, so each half reads standalone
# waterfowl is measured as an absolute level in BOTH models - there is no anomaly version of
# it - so it belongs in row (a) and is simply absent from row (b)
fig3 <- panel_plot(both[kind == "level"], "a. Absolute conditions") /
        panel_plot(both[kind == "shock" & variable != "Anseriformes"],
                   "b. Anomalies (departure from the 90-day baseline)")
ggsave_spark(file.path(figures_main_folder, "figure3_levels_vs_shocks.png"),
             fig3, width = 13, height = 8.5)

# Figure 4 lives in d_17 now - it's the bird cumulative curves over 0-28 vs 0-56 d, which needs
# its own fits (per-variable lag windows, predators as a comparator). Nothing here builds it.
saveRDS(list(fig3_data = both),
        paste0(objects_folder, "figures_3_4_data.RDS"))
cat("wrote figure3_levels_vs_shocks.png\n")
