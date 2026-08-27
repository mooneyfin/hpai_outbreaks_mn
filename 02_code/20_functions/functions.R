############ DATA PREPARATION FUNCTIONS ############

############ ANALYSIS FUNCTIONS ############

############ FIGURE PALETTES AND THEME HELPERS ############

#1a.Sequential rate palette (low → high, 4 levels)
spark_sequential_rate_colors <- c("#F2EBCF", "#C4AEDB", "#6B4A9A", "#32135C")

#1b.Diverging rate palette (anchored at 1.0)
spark_diverging_rate_colors <- c(low = "#2D6FA3", mid = "#F2EBCF", high = "#C9405A")

#1c.Sex palette (Female first, alphabetical)
spark_sex_colors <- c("Female" = "#3F6961", "Male" = "#EAB63E")

#1d.Age-group palette (18 levels — use full palette even when plotting a subset)
spark_age_group_colors <- c(
  "#770042", "#8E0052", "#AA0265", "#CB0A7A", "#E1218A", "#EB3A94",
  "#EF549E", "#F06EA9", "#F495BF", "#AA9FCE", "#61CDF6", "#00B9F0",
  "#00ACE8", "#00A2DD", "#0094CB", "#007CAA", "#005A7E", "#003A57"
)

#2a.Standard theme for non-map plots (line, bar, forest)
theme_spark <- function(base_size = 11) {
  ggplot2::theme_classic(base_family = "Avenir", base_size = base_size) +
    ggplot2::theme(
      axis.text.x      = ggplot2::element_text(angle = 45, hjust = 1),
      axis.line        = ggplot2::element_line(linewidth = 0.25 / ggplot2::.pt),
      axis.ticks       = ggplot2::element_line(linewidth = 0.25 / ggplot2::.pt),
      plot.title       = ggplot2::element_text(family = "Avenir Heavy", hjust = 0),
      strip.background = ggplot2::element_blank(),
      strip.text       = ggplot2::element_text(family = "Avenir Heavy", hjust = 0, size = ggplot2::rel(1)),
      legend.position  = "bottom",
      legend.key.size  = ggplot2::unit(8, "pt")
    )
}

#2b.Theme for maps (transform sf objects to EPSG:5070 before plotting US maps)
# ggsave through ragg, otherwise the Avenir faces silently fall back to the default.
ggsave_spark <- function(file, plot, width, height, dpi = 300, bg = "white") {
  # bg has to be explicit: theme_void leaves plot.background blank, so map and flowchart
  # exports come out transparent and journals render them on black
  dev <- if (requireNamespace("ragg", quietly = TRUE)) ragg::agg_png else NULL
  ggplot2::ggsave(file, plot, width = width, height = height, dpi = dpi, device = dev, bg = bg)
}

theme_spark_map <- function() {
  ggplot2::theme_void(base_family = "Avenir") +
    ggplot2::theme(
      plot.title          = ggplot2::element_text(family = "Avenir Heavy", hjust = 0,
                                                  size = 12, margin = ggplot2::margin(b = 8)),
      plot.subtitle       = ggplot2::element_text(family = "Avenir", hjust = 0,
                                                  size = 11, margin = ggplot2::margin(b = 8),
                                                  color = "grey30"),
      plot.margin         = ggplot2::margin(t = 12, r = 5, b = 5, l = 5),
      strip.text          = ggplot2::element_text(family = "Avenir Heavy", hjust = 0, size = 11),
      legend.position     = "bottom",
      legend.justification = "centre",
      legend.margin       = ggplot2::margin(l = 8),
      legend.title        = ggplot2::element_text(family = "Avenir", size = 7),
      legend.text         = ggplot2::element_text(family = "Avenir", size = 7)
    )
}
