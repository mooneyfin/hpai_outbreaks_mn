############ FIGURE 2 — CUMULATIVE RATE RATIOS, BOTH PRIMARY MODELS ############
# The version of this figure that was submitted had no source script, only a PNG, and it had
# dropped Anseriformes from the anomaly panel - which is where the paper's strongest result
# lives (Minnesota, lag 22-28 d, RR 2.53). Both panels now carry all five exposures.
#
# Waterfowl abundance enters BOTH models in absolute form; only the four meteorological terms
# are re-expressed as 90-day shocks. That is worth a line in the caption,
# otherwise panel b reads as though the bird term were anomalised too.
#
# Nothing is refitted: the cumulative tables come straight off the cached fits.

project.folder = paste0(print(here::here()), '/')
src <- readLines(paste0(project.folder, 'create_folder_structure.R'))
eval(parse(text = paste(src[!grepl("^\\s*rm\\(list", src)], collapse = "\n")))
source(paste0(functions.folder, 'script_initiate.R'))
suppressMessages({library(viridisLite); library(patchwork)})
source(paste0(functions.folder, 'inla_dlnm_helpers.R'))   # exposure set and panel labels

PANEL_ORD <- INLA_PANEL_ORD
FITKEY    <- setNames(c("Minnesota", "Other flyway", "Pooled"), INLA_PANEL_ORD)
# exposures come from INLA_PRIMARY_* in inla_dlnm_helpers.R - never redefined locally
EXPO_ORD  <- unname(INLA_PRIMARY_LAB[c(INLA_PRIMARY_MET, INLA_BIRD_VAR)])
MODELS    <- c(level = "a. Absolute conditions",
               shock = "b. 90-day shocks (departure from the preceding 90 days)")
# the single-panel main figure gets a plain title in place of the panel legend
OK  <- setNames(viridisLite::viridis(3, begin = 0.05, end = 0.78), PANEL_ORD)
SHP <- setNames(c(16, 17, 15), PANEL_ORD)
lab_num <- scales::label_number(drop0trailing = TRUE)

fits <- readRDS(paste0(objects_folder, "primary_shock_level_fits.RDS"))

# "1.24 (0.82, 1.87)" -> three numbers. Fails loudly rather than silently returning NA, since a
# quietly dropped row is exactly how Anseriformes went missing last time.
parse_est <- function(x) {
  g <- regmatches(x, regexec("^([0-9.]+) \\(([0-9.]+), ([0-9.]+)\\)$", x))
  stopifnot(all(lengths(g) == 4))
  as.data.table(setNames(lapply(2:4, function(i) as.numeric(vapply(g, `[`, "", i))),
                         c("RR", "low", "high")))
}

dat <- rbindlist(lapply(names(MODELS), function(mk) rbindlist(lapply(PANEL_ORD, function(pn) {
  k <- paste0(FITKEY[[pn]], "|", mk)
  stopifnot(k %in% names(fits))
  cu <- as.data.table(fits[[k]]$cumulative)
  wcols <- setdiff(names(cu), "Predictor")
  rbindlist(lapply(wcols, function(w)
    cbind(data.table(model = MODELS[[mk]], panel = pn, n = fits[[k]]$n,
                     variable = cu$Predictor, window = w),
          parse_est(cu[[w]]))))
}))))

# every exposure must survive into every panel, or the figure is lying by omission
stopifnot(!anyNA(dat$RR), all(EXPO_ORD %in% dat$variable))
chk <- dat[, .N, by = .(model, panel, window)]
stopifnot(all(chk$N == length(EXPO_ORD)))
cat(sprintf("rows %d | models %d | panels %d | exposures %d\n", nrow(dat),
            uniqueN(dat$model), uniqueN(dat$panel), uniqueN(dat$variable)))

WIN_ORD <- unique(dat$window)
dat[, `:=`(panel    = factor(panel, levels = PANEL_ORD),
           window   = factor(window, levels = WIN_ORD),
           variable = factor(variable, levels = rev(EXPO_ORD)),
           model    = factor(model, levels = unname(MODELS)))]

# The panel label rides in the LEGEND TITLE slot rather than plot.title. That is what puts
# "a. Absolute conditions" on the same line as the three keys instead of stacking it above them,
# and it means each panel carries its own legend without the two ever drifting apart.
one_panel <- function(mk, keep_xlab) {
  ggplot(dat[model == MODELS[[mk]]], aes(RR, variable, colour = panel, shape = panel)) +
    geom_vline(xintercept = 1, colour = "grey45", linewidth = 0.3) +
    geom_errorbarh(aes(xmin = low, xmax = high), height = 0, linewidth = 0.55,
                   position = position_dodge(width = 0.66, reverse = TRUE)) +
    geom_point(size = 2.1, position = position_dodge(width = 0.66, reverse = TRUE)) +
    facet_wrap(~ window, nrow = 1) +
    scale_colour_manual(values = OK) + scale_shape_manual(values = SHP) +
    # log axis: a ratio of 2 and a ratio of 0.5 are the same size of effect, and on a linear
    # axis they are not. coord_cartesian, never scale_x limits, which DELETES out-of-range rows.
    scale_x_continuous(trans = "log", breaks = c(0.5, 1, 2, 3), labels = lab_num) +
    coord_cartesian(xlim = c(0.5, 4.2)) +
    labs(x = if (keep_xlab) "Rate ratio per +0.5 standard deviation increase in exposure" else NULL,
         y = NULL,
         colour = if (uniqueN(dat$panel) == 1) "Minnesota, 90-day shocks" else MODELS[[mk]],
         shape  = if (uniqueN(dat$panel) == 1) "Minnesota, 90-day shocks" else MODELS[[mk]]) +
    guides(colour = guide_legend(title.position = "left"),
           shape  = guide_legend(title.position = "left")) +
    theme_classic(base_family = "Avenir") +
    theme(legend.position = "top", legend.justification = "left",
          legend.text = element_text(size = 12),
          legend.title = element_text(family = "Avenir Heavy", size = 13, vjust = 0.5),
          legend.margin = margin(b = 2),
          strip.text = element_text(family = "Avenir Heavy", size = 11),
          strip.background = element_blank(),
          axis.text = element_text(size = 10, colour = "black"),
          axis.title = element_text(size = 11),
          panel.border = element_rect(fill = NA, colour = "grey20", linewidth = 0.4))
}

# Figure 2 (main text): Minnesota, 90-day shock - the primary model on its own.
# Figure S2: the full grid, both parameterisations, all three panels, for the supplement.
dat_all <- copy(dat)
dat <- dat_all[panel == "Minnesota" & model == MODELS[["shock"]]]
# one panel, one parameterisation: nothing for a legend to distinguish, so it goes
f2 <- one_panel("shock", TRUE) + theme(legend.position = "none")
ggsave_spark(file.path(figures_main_folder, "figure2_primary_forest.png"), f2, width = 13, height = 4.2)
dat <- dat_all
fS2 <- one_panel("level", FALSE) / one_panel("shock", TRUE) + plot_layout(heights = c(1, 1))
ggsave_spark(file.path(figures_supplement_folder, "figureS2_primary_forest_all_panels.png"), fS2, width = 13, height = 9)
saveRDS(dat_all, paste0(objects_folder, "figure2_primary_forest.RDS"))
fwrite(dat_all, paste0(objects_folder, "figure2_primary_forest.csv"))
cat("wrote figure2_primary_forest.png (Minnesota, shock) and figureS2_primary_forest_all_panels.png\n")
print(dcast(dat_all[variable == "Anseriformes"], model + panel ~ window, value.var = "RR"))
