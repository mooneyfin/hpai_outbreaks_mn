############ LAG-RESPONSE CURVES ############
# Cumulative OR accumulated from day 0 out to each lag, Minnesota only, every variable.
# Same construction as the bird figure, so the two read the same way.
# (a) absolute conditions, (b) anomalies.

rm(list = ls())
project.folder = paste0(print(here::here()), '/')
source(paste0(project.folder, 'create_folder_structure.R'))
source(paste0(functions.folder, 'script_initiate.R'))
suppressMessages({library(INLA); library(Matrix); library(patchwork); library(viridisLite)})
source(paste0(functions.folder, 'inla_dlnm_helpers.R'))

fits    <- readRDS(paste0(objects_folder, "primary_shock_level_fits.RDS"))
lab_num <- scales::label_number(drop0trailing = TRUE)
LINE    <- viridisLite::viridis(1, begin = 0.05)

# 1. two readings of the same fit: accumulated from day 0 out to each lag, or the effect of
#    that single lag day on its own
er_curve <- function(r, v, kind) {
  rbindlist(lapply(0:r$obj$max_lag, function(a) {
    lags <- if (kind == "cumulative") 0:a else a
    e <- met_lag_effect(r$fit, r$obj, v, lags)
    data.table(variable = unname(r$labels[v]), x = a,
               OR = e[["OR"]], low = e[["low"]], high = e[["high"]])
  }))
}

build <- function(key, kind, drop_bird = FALSE) {
  r <- fits[[key]]
  vs <- r$obj$met_vars
  if (drop_bird) vs <- setdiff(vs, "anseriformes")
  d <- rbindlist(lapply(vs, er_curve, r = r, kind = kind))
  d[, variable := sub("^Surface runoff$", "Runoff", variable)]
  d[, variable := sub(" shock$", "", variable)]
  d
}

ORD <- c("Temperature", "Precipitation", "Soil moisture", "Wind speed", "Anseriformes")
YCAP <- 6
prep <- function(kind) {
  a <- build("Minnesota|level", kind)
  b <- build("Minnesota|shock", kind, drop_bird = TRUE)   # no anomaly version of waterfowl
  for (d in list(a, b)) {
    d[, variable := factor(variable, levels = ORD)]
    # clip so one wide interval doesn't flatten every other panel
    d[, `:=`(low = pmax(low, 1 / YCAP), high = pmin(high, YCAP))]
  }
  list(abs = a, ano = b)
}

er_plot <- function(d, ttl, xlab) {
  ggplot(d, aes(x, OR)) +
    geom_hline(yintercept = 1, colour = "grey45", linewidth = 0.3) +
    geom_ribbon(aes(ymin = low, ymax = high), fill = LINE, alpha = 0.16) +
    geom_line(colour = LINE, linewidth = 0.95) +
    facet_wrap(~ variable, nrow = 1, drop = TRUE) +
    scale_y_continuous(trans = "log", breaks = c(0.25, 0.5, 1, 2, 4), labels = lab_num) +
    scale_x_continuous(breaks = seq(0, 28, 7)) +
    labs(x = xlab, y = "Odds ratio", title = ttl) +
    theme_classic(base_family = "Avenir") +
    theme(plot.title = element_text(family = "Avenir Heavy", size = 14,
                                    margin = margin(b = 4)),
          plot.title.position = "plot",
          strip.text = element_text(family = "Avenir Heavy", size = 11),
          strip.background = element_blank(),
          axis.title = element_text(size = 12),
          axis.text = element_text(size = 10, colour = "black"),
          panel.border = element_rect(fill = NA, colour = "grey20", linewidth = 0.4))
}

out <- list()
for (kd in c("cumulative", "lagspecific")) {
  p <- prep(kd)
  fig <- er_plot(p$abs, "a. Absolute conditions", "Lag (days)") /
         er_plot(p$ano, "b. Anomalies (departure from the 90-day baseline)", "Lag (days)")
  f <- sprintf("figure5_lag_response_%s.png", kd)
  ggsave_spark(file.path(figures_main_folder, f), fig, width = 13, height = 7.5)
  out[[kd]] <- p
  cat("wrote ", f, "\n", sep = "")
}
saveRDS(out, paste0(objects_folder, "figure_lag_response_data.RDS"))
