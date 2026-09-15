############ EPIDEMIC CURVE ############
# Weekly spillover events, Minnesota stacked against the other northern Mississippi Flyway
# states. Meteorological season boundaries marked (Mar 1 / Jun 1 / Sep 1 / Dec 1), which is
# how `season` is defined everywhere else in the pipeline.
#
# Flyway membership follows the USFWS administrative flyways: the seven states in the panel
# (MN, WI, IA, MI, IL, IN, OH) are the northern Mississippi Flyway.

rm(list = ls())
project.folder = paste0(print(here::here()), '/')
source(paste0(project.folder, 'create_folder_structure.R'))
source(paste0(functions.folder, 'script_initiate.R'))
suppressMessages(library(viridisLite))

panel <- setDT(readRDS(paste0(objects_folder, multistate_ts_daily_rds)))
scr   <- readRDS(paste0(objects_folder, multistate_independence_screening_rds))

# same case definition as the models: drop cross-farm CSLT/NOSEQ events
count_flagged <- function(s, ids) {
  if (is.na(s)) return(0L)
  as.integer(sum(as.numeric(trimws(unlist(strsplit(s, ",")))) %in% ids))
}
panel[, nf := vapply(outbreak_id, count_flagged, integer(1), ids = scr$cross_farm_ids)]
panel[, oc := pmax(outbreak_count - nf, 0L)]

GRP <- c("Minnesota", "Other northern Mississippi Flyway states")
cases <- panel[oc > 0]
cases[, grp := factor(fifelse(state == "Minnesota", GRP[1], GRP[2]), levels = GRP)]
cases[, week := as.integer(strftime(date, "%V"))]
wk <- cases[, .(n = sum(oc)), by = .(week, grp)]
cat(sprintf("%d events total | Minnesota %d | other flyway %d\n",
            sum(wk$n), sum(wk[grp == GRP[1]]$n), sum(wk[grp == GRP[2]]$n)))

# meteorological season starts as ISO weeks
seas <- data.table(week = c(9, 22, 35, 48),
                   label = c("Spring", "Summer", "Fall", "Winter"))
# two-line ticks: week number always, month name on the first tick falling in each month
brk     <- seq(2, 52, 2)
brk_mon <- strftime(as.Date("2022-01-01") + (brk - 1) * 7, "%b")
brk_lab <- paste0(brk, "\n", ifelse(c(TRUE, diff(match(brk_mon, month.abb)) != 0), brk_mon, ""))

COL   <- setNames(viridisLite::viridis(2, begin = 0.05, end = 0.62), GRP)
ymax  <- wk[, sum(n), by = week][, max(V1)]

fig <- ggplot(wk, aes(week, n, fill = grp)) +
  geom_col(colour = "black", linewidth = 0.25, width = 0.9) +
  geom_vline(data = seas, aes(xintercept = week), colour = "red",
             linetype = "dashed", linewidth = 0.7) +
  geom_text(data = seas, aes(x = week - 0.6, y = ymax * 0.98, label = label),
            inherit.aes = FALSE, colour = "red", family = "Avenir Heavy", size = 3.6,
            angle = 90, hjust = 1, vjust = 0.5) +
  scale_fill_manual(values = COL) +
  scale_x_continuous(breaks = brk, limits = c(0.5, 52.5), expand = c(0, 0), labels = brk_lab) +
  scale_y_continuous(breaks = scales::pretty_breaks(8), expand = expansion(mult = c(0, 0.05)),
                     labels = scales::label_number(drop0trailing = TRUE)) +
  labs(x = "Week", y = "Number of spillover events", fill = NULL) +
  theme_classic(base_family = "Avenir") +
  theme(panel.border = element_rect(fill = NA, colour = "grey20", linewidth = 0.5),
        legend.position = "top", legend.justification = "left",
        legend.text = element_text(size = 12),
        legend.key.size = unit(1.3, "lines"),
        axis.title = element_text(size = 12),
        axis.text  = element_text(size = 10, colour = "black"))

ggsave_spark(file.path(figures_intermediate_folder, "epidemic_curve.png"),
             fig, width = 12, height = 5.5)
saveRDS(wk, paste0(objects_folder, "figure_epidemic_curve_data.RDS"))
cat("wrote intermediate/epidemic_curve.png\n")
