############ FIGURE 4 — METEOROLOGY ONLY MATTERS IN THE WAKE OF WATERFOWL ############
# The meteorological effect at lag 0-7 and 8-14 days, shown separately for cell-days with high and low
# Anseriformes abundance three weeks earlier. This is the pre-stated interaction model (c_50):
# strata lag basis, continuous modifier, primary prior, both parameterisations, all panels.
#
# Both pre-stated windows are drawn, 0-7 and 8-14 days. The hypothesis named "the two weeks
# before", and showing only the window that produced the result would be choosing it after the
# fact. Anseriformes is not a row here: it is the modifier, not an exposure being modified.
#
# Nothing is refitted: the estimates are read off the cached c_50 output. The parse is strict so
# a missing exposure fails loudly rather than vanishing from the figure.

project.folder = paste0(print(here::here()), '/')
src <- readLines(paste0(project.folder, 'create_folder_structure.R'))
eval(parse(text = paste(src[!grepl("^\\s*rm\\(list", src)], collapse = "\n")))
source(paste0(functions.folder, 'script_initiate.R'))
source(paste0(functions.folder, 'inla_dlnm_helpers.R'))
suppressMessages(library(patchwork))

PANEL_ORD <- c("Minnesota", "Other flyway states", "Pooled")
EXPO_ORD  <- unname(INLA_PRIMARY_LAB[INLA_PRIMARY_MET])
MODELS    <- c(`Absolute conditions` = "a. Absolute conditions",
               `90-day shock`        = "b. 90-day shocks (departure from the preceding 90 days)")
GRP <- c(`High waterfowl` = "High waterfowl (+1 SD, lag 21-28 d)",
         `Low waterfowl`  = "Low waterfowl (-1 SD, lag 21-28 d)")
COL <- setNames(c("#D55E00", "#0072B2"), unname(GRP))
SHP <- setNames(c(16, 17), unname(GRP))
lab_num <- scales::label_number(drop0trailing = TRUE)

ip <- setDT(readRDS(paste0(objects_folder, "si_interaction_primary.RDS")))
ip[, Kind := fifelse(Kind == "90-day anomaly", "90-day shock", Kind)]

# "1.64 (1.12, 2.44)" -> three numbers, or fail
parse_est <- function(x) {
  g <- regmatches(x, regexec("^([0-9.]+) \\(([0-9.]+), ([0-9.]+)\\)$", x))
  stopifnot(all(lengths(g) == 4))
  as.data.table(setNames(lapply(2:4, function(i) as.numeric(vapply(g, `[`, "", i))),
                         c("RR", "low", "high")))
}
dat <- rbindlist(lapply(names(GRP), function(g)
  cbind(ip[, .(Kind, Panel, Exposure, Window, n)], group = GRP[[g]], parse_est(ip[[g]]))))
dat <- dat[Exposure %in% EXPO_ORD]
stopifnot(!anyNA(dat$RR), all(EXPO_ORD %in% dat$Exposure), all(PANEL_ORD %in% dat$Panel))
dat[, `:=`(Panel = factor(Panel, levels = PANEL_ORD),
           Window = factor(paste("Lag", Window), levels = c("Lag 0-7 days", "Lag 8-14 days")),
           Exposure = factor(Exposure, levels = rev(EXPO_ORD)),
           group = factor(group, levels = unname(GRP)))]
fwrite(dat, paste0(objects_folder, "figure4_interaction.csv"))

one_block <- function(kind, keep_xlab) {
  ggplot(dat[Kind == kind], aes(RR, Exposure, colour = group, shape = group)) +
    geom_vline(xintercept = 1, colour = "grey45", linewidth = 0.3) +
    geom_errorbar(aes(xmin = low, xmax = high), width = 0, linewidth = 0.55,
                  orientation = "y", position = position_dodge(width = 0.6, reverse = TRUE)) +
    geom_point(size = 2.2, position = position_dodge(width = 0.6, reverse = TRUE)) +
    facet_grid(Window ~ Panel) +
    scale_colour_manual(values = COL) + scale_shape_manual(values = SHP) +
    # log axis so a halving and a doubling are the same distance; coord_cartesian clips
    # rather than deleting any estimate whose interval runs past the edge
    scale_x_continuous(trans = "log", breaks = c(0.25, 0.5, 1, 2, 3), labels = lab_num) +
    coord_cartesian(xlim = c(0.18, 3.6)) +   # 0.36 (0.20, 0.67) sits at the left edge
    labs(x = if (keep_xlab) "Rate ratio per +0.5 SD increase in exposure" else NULL,
         y = NULL, colour = MODELS[[kind]], shape = MODELS[[kind]]) +
    guides(colour = guide_legend(title.position = "left"), shape = guide_legend(title.position = "left")) +
    theme_classic(base_family = "Avenir") +
    theme(legend.position = "top", legend.justification = "left",
          legend.text = element_text(size = 11),
          legend.title = element_text(family = "Avenir Heavy", size = 13, vjust = 0.5),
          legend.margin = margin(b = 2),
          strip.text = element_text(family = "Avenir Heavy", size = 11),
          strip.background = element_blank(),
          axis.text = element_text(size = 10, colour = "black"),
          axis.title = element_text(size = 11),
          panel.border = element_rect(fill = NA, colour = "grey20", linewidth = 0.4))
}
fig <- one_block("Absolute conditions", FALSE) / one_block("90-day shock", TRUE)
ggsave_spark(file.path(figures_main_folder, "figure4_interaction.png"), fig, width = 12, height = 13)
cat("wrote figure4_interaction.png\n")
print(dcast(dat[Exposure == "Temperature"], Kind + Panel + Window ~ group, value.var = "RR"))
