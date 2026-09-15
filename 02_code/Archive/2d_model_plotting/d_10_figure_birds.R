############ FIGURE 4 — BIRD-ONLY MODEL ############
# Anseriformes and predators in their own model, no meteorology, over a long lag window.
#
# Cumulative only, deliberately. The lag-SPECIFIC profile is not identified here: on
# identical data it rises with lag under a 0-28 d window and falls under 0-56 d, and WAIC
# separates those by 1.8. The cumulative is stable across the same variants (mean per-day OR
# 1.049 vs 1.050), so it is the part the data actually supports.
#
# Two things about the panel choice. This model uses `case_crossover_df_timestrat_month_post7`
# rather than the shock panel: with no shock terms there is no 90-day warm-up to satisfy, so
# it keeps all 175 cases instead of 167. And it runs at 0-56 days because the whole point is
# the long-lag migration signal - the 28-day cap exists for the meteorological terms.
#
# Pooled is deliberately absent: Minnesota and the other flyway states are the two
# independent samples, and pooling them here would just average a result the figure is
# meant to show replicating.

rm(list = ls())
project.folder = paste0(print(here::here()), '/')
source(paste0(project.folder, 'create_folder_structure.R'))
source(paste0(functions.folder, 'script_initiate.R'))
suppressMessages({library(INLA); library(Matrix); library(viridisLite)})
source(paste0(functions.folder, 'inla_dlnm_helpers.R'))

PREC    <- 4                # N(0, 0.5^2), same prior as the primary
MAXLAG  <- 56
BIRDS   <- c("anseriformes", "predators")
BLAB    <- c(anseriformes = "Anseriformes", predators = "Predator birds")
PANELS  <- c("Minnesota", "Other flyway states")
COL     <- setNames(viridisLite::viridis(2, begin = 0.05, end = 0.62), PANELS)
lab_num <- scales::label_number(drop0trailing = TRUE)

D <- setDT(readRDS(paste0(objects_folder, "case_crossover_df_timestrat_month_post7.RDS")))
if (anyDuplicated(names(D))) D <- D[, which(!duplicated(names(D))), with = FALSE]
sub <- list(Minnesota = D[state == "Minnesota"],
            `Other flyway states` = D[state != "Minnesota"])

av <- setNames(rep(list(list(fun = "lin")), length(BIRDS)), BIRDS)
res <- rbindlist(lapply(PANELS, function(pn) {
  d <- sub[[pn]]
  o <- assemble_cc_dlnm(d, met_vars = BIRDS, weekly_vars = character(0), max_lag = MAXLAG,
                        argvar_by_var = av, arglag = list(fun = "ns", df = 2))
  f <- suppressWarnings(fit_cc_inla_dlnm(o, fixed_prec = PREC))
  cat(sprintf("%-13s n = %3d | WAIC %.1f\n", pn, sum(d$outbreak_binary), f$waic$waic))
  rbindlist(lapply(BIRDS, function(v) rbindlist(lapply(0:MAXLAG, function(a) {
    cu <- met_lag_effect(f, o, v, 0:a)        # accumulated from day 0 to that day
    data.table(panel = pn, taxon = BLAB[[v]], lag = a,
               OR = cu[["OR"]], low = cu[["low"]], high = cu[["high"]])
  }))))
}))
res[, panel := factor(panel, levels = PANELS)]
res[, taxon := factor(taxon, levels = unname(BLAB))]

fig <- ggplot(res, aes(lag, OR, colour = panel, fill = panel)) +
  geom_hline(yintercept = 1, colour = "grey45", linewidth = 0.3) +
  geom_ribbon(aes(ymin = low, ymax = high), alpha = 0.13, colour = NA) +
  geom_line(aes(linetype = panel), linewidth = 0.95) +
  facet_wrap(~ taxon, nrow = 1) +
  scale_colour_manual(values = COL) + scale_fill_manual(values = COL) +
  scale_linetype_manual(values = setNames(c("solid", "dashed"), PANELS)) +
  scale_x_continuous(breaks = seq(0, MAXLAG, 14)) +
  # shared log y: anseriformes runs 1.0-11.4 and predators 0.3-4.1, so on a linear shared
  # axis predators flattens to a line. on log both are legible and OR = 1 lands at the same
  # height in both panels, which is the comparison the figure is for.
  scale_y_continuous(trans = "log", breaks = c(0.25, 0.5, 1, 2, 4, 8), labels = lab_num) +
  labs(x = "Lag (days)",
       y = "Odds ratio per +0.5 standard deviation increase in abundance",
       colour = NULL, fill = NULL, linetype = NULL) +
  theme_classic(base_family = "Avenir") +
  # legend tucked into the top-left of the Predators panel, which is empty
  theme(legend.position = c(0.99, 0.97), legend.justification = c(1, 1),
        legend.text = element_text(size = 13),
        legend.key.size = unit(1.4, "lines"),
        legend.background = element_rect(fill = "white", colour = NA),
        strip.text = element_text(family = "Avenir Heavy", size = 13),
        strip.background = element_blank(),
        axis.title = element_text(size = 12),
        axis.text = element_text(size = 10, colour = "black"),
        panel.border = element_rect(fill = NA, colour = "grey20", linewidth = 0.4))

ggsave_spark(file.path(figures_main_folder, "figure4_birds_lag.png"), fig,
             width = 10, height = 5.5)
saveRDS(res, paste0(objects_folder, "figure4_bird_model_data.RDS"))

# the numbers behind the long-lag windows, for the text
cat("\ncumulative OR at selected lags:\n")
print(dcast(res[lag %in% c(7, 14, 28, 42, 56)],
            taxon + lag ~ panel,
            value.var = "OR", fun.aggregate = function(x) sprintf("%.2f", x[1])))
cat("\nwrote figure4_birds_lag.png\n")
