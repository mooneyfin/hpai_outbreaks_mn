############ FIGURE 2 — STUDY SETTING ############
# (a) case-cell map across the seven states
# (b) epidemic curve, Minnesota stacked against the other northern Mississippi Flyway states
#
# The map is built in Python (geopandas + contextily, d_02_case_zone_map.py) because that is
# where the basemap tooling lives, so it comes in here as a raster rather than a ggplot.
# Run d_01 and d_02 before this.

rm(list = ls())
project.folder = paste0(print(here::here()), '/')
source(paste0(project.folder, 'create_folder_structure.R'))
source(paste0(functions.folder, 'script_initiate.R'))
suppressMessages({library(patchwork); library(magick); library(viridisLite)})

wk  <- setDT(readRDS(paste0(objects_folder, "figure_epidemic_curve_data.RDS")))
GRP <- levels(wk$grp)

# 1b. Weekly Anseriformes abundance over the case cells, to overlay on the epidemic curve. The
#     point is that migration LEADS spillover: waterfowl peak week 11 in Minnesota and week 9
#     elsewhere, while the outbreak curve peaks weeks 13-17.
#     Same ISO-week definition d_01 used for the bars, or the two would not line up.
dailyp <- setDT(readRDS(paste0(objects_folder, multistate_ts_daily_rds)))
dailyp[, date := as.Date(date)]
ccp <- setDT(readRDS(paste0(objects_folder, "case_crossover_df_shock_post7.RDS")))
if (anyDuplicated(names(ccp))) ccp <- ccp[, which(!duplicated(names(ccp))), with = FALSE]
case_cells <- unique(ccp[outbreak_binary == 1]$zone_id)

abund <- dailyp[zone_id %in% case_cells & date >= as.Date("2022-01-01"),
                .(ans = mean(anseriformes, na.rm = TRUE)),
                by = .(week = as.integer(strftime(date, "%V")),
                       grp = factor(fifelse(state == "Minnesota", GRP[1], GRP[2]), levels = GRP))]
# abundance runs to ~93, events to 20, so the secondary axis is events = abundance / SCALE
SCALE <- 5
COL <- setNames(viridisLite::viridis(2, begin = 0.05, end = 0.62), GRP)

seas <- data.table(week = c(9, 22, 35, 48), label = c("Spring", "Summer", "Fall", "Winter"))
brk     <- seq(2, 52, 2)
brk_mon <- strftime(as.Date("2022-01-01") + (brk - 1) * 7, "%b")
brk_lab <- paste0(brk, "\n", ifelse(c(TRUE, diff(match(brk_mon, month.abb)) != 0), brk_mon, ""))
ymax <- wk[, sum(n), by = week][, max(V1)]

pa <- ggplot(wk, aes(week, n, fill = grp)) +
  geom_col(colour = "grey35", linewidth = 0.2, width = 0.9, alpha = 0.45) +
  geom_vline(data = seas, aes(xintercept = week), colour = "grey45",
             linetype = "dashed", linewidth = 0.5) +
  geom_text(data = seas, aes(x = week - 0.6, y = ymax * 0.98, label = label),
            inherit.aes = FALSE, colour = "grey35", family = "Avenir Heavy", size = 4.6,
            angle = 90, hjust = 1, vjust = 0.5) +
  geom_line(data = abund, aes(week, ans / SCALE, group = grp),
            inherit.aes = FALSE, colour = "white", linewidth = 2.6, lineend = "round") +
  geom_line(data = abund, aes(week, ans / SCALE, colour = grp),
            inherit.aes = FALSE, linewidth = 1.4, lineend = "round") +
  scale_fill_manual(values = COL) +
  scale_colour_manual(values = COL, guide = "none") +
  scale_x_continuous(breaks = brk, limits = c(0.5, 52.5), expand = c(0, 0), labels = brk_lab) +
  scale_y_continuous(breaks = scales::pretty_breaks(8), expand = expansion(mult = c(0, 0.05)),
                     labels = scales::label_number(drop0trailing = TRUE),
                     sec.axis = sec_axis(~ . * SCALE, name = "Anseriformes abundance (lines)",
                                         breaks = scales::pretty_breaks(6))) +
  labs(x = "Week", y = "Number of spillover events", fill = "b. ") +
  theme_classic(base_family = "Avenir") +
  theme(panel.border = element_rect(fill = NA, colour = "grey20", linewidth = 0.5),
        legend.position = "top", legend.justification = "left",
        legend.title = element_text(family = "Avenir Heavy", size = 16),
        legend.text = element_text(size = 15), legend.key.size = unit(1.5, "lines"),
        axis.title = element_text(size = 15),
        axis.text = element_text(size = 13, colour = "black"),
        axis.title.y.right = element_text(size = 15, margin = margin(l = 10)))

mp  <- magick::image_read(file.path(figures_intermediate_folder, "case_zone_map.png"))
asp <- with(magick::image_info(mp), height / width)

# rasterGrob at full npc so the map fills its panel edge to edge. that only looks right if the
# panel has the map's own aspect, which is what the heights below do - otherwise it stretches.
pmap <- ggplot() +
  annotation_custom(grid::rasterGrob(as.raster(mp), interpolate = TRUE,
                                     width = unit(1, "npc"), height = unit(1, "npc"))) +
  labs(title = "a. ") +
  theme_void(base_family = "Avenir") +
  theme(plot.title = element_text(family = "Avenir Heavy", size = 16, hjust = 0),
        plot.margin = margin(0, 0, 6, 0))

# same width for both panels: give the map exactly W * aspect inches of height and let the
# curve take the rest, so neither is letterboxed or squeezed
W       <- 11
H_CURVE <- 4.6
H_MAP   <- W * asp
fig <- pmap / pa + plot_layout(heights = c(H_MAP, H_CURVE))
ggsave_spark(file.path(figures_main_folder, "figure1_study_setting.png"),
             fig, width = W, height = H_MAP + H_CURVE)
cat("wrote figure1_study_setting.png\n")
