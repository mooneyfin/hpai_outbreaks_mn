# d_08_consort_flowchart.R
# Sample-size cascade for the daily case-crossover study of HPAI spillover into poultry across
# the seven northern Mississippi Flyway states.
#
# Starts at the full flyway panel, screens out farm-to-farm spread, then branches into Minnesota
# and the other six states because those are the two panels every model is fitted on. The
# branches rejoin at the pooled analytic sample.
#
# Every number is counted from the panels themselves rather than hardcoded, so it can't drift
# away from the models the way the old cascade objects did.

# 0a. Root, folders, packages
rm(list = ls())
project.folder = paste0(print(here::here()), '/')
source(paste0(project.folder, 'create_folder_structure.R'))
source(paste0(functions.folder, 'script_initiate.R'))

fmt <- function(x) formatC(x, big.mark = ",", format = "d")

# 1a. Full flyway panel, before anything is excluded
panel <- setDT(readRDS(paste0(objects_folder, multistate_ts_daily_rds)))
scr   <- readRDS(paste0(objects_folder, multistate_independence_screening_rds))
n_states    <- uniqueN(panel$state)
n_celldays  <- nrow(panel)
n_gridcells <- uniqueN(panel$zone_id)
date_from   <- format(min(panel$date), "%d %b %Y")
date_to     <- format(max(panel$date), "%d %b %Y")

n_raw_events   <- sum(panel$outbreak_count, na.rm = TRUE)
n_raw_celldays <- sum(panel$outbreak_count > 0, na.rm = TRUE)
n_raw_cells    <- uniqueN(panel[outbreak_count > 0]$zone_id)

# 1b. Independence screen: events flagged by sequence and epidemiological investigation as
#     farm-to-farm spread aren't independent introductions from the environment
count_flagged <- function(s, ids) {
  if (is.na(s)) return(0L)
  as.integer(sum(as.numeric(trimws(unlist(strsplit(s, ",")))) %in% ids))
}
panel[, nf := vapply(outbreak_id, count_flagged, integer(1), ids = scr$cross_farm_ids)]
n_flagged <- sum(panel$nf)
panel[, oc := pmax(outbreak_count - nf, 0L)]
n_elig_events   <- sum(panel$oc)
n_elig_celldays <- sum(panel$oc > 0)
n_elig_cells    <- uniqueN(panel[oc > 0]$zone_id)

elig <- panel[oc > 0][, .(events = sum(oc), cells = uniqueN(zone_id)),
                      by = .(g = fifelse(state == "Minnesota", "mn", "other"))]
E <- function(grp, f) elig[g == grp][[f]]

# 1c. The two case-crossover panels. `pre` is the design before the 90-day anomaly needs a
#     baseline; `cc` is what the primary models are actually fitted to.
load_cc <- function(f) {
  d <- setDT(readRDS(paste0(objects_folder, f)))
  if (anyDuplicated(names(d))) d <- d[, which(!duplicated(names(d))), with = FALSE]
  d[, g := fifelse(state == "Minnesota", "mn", "other")][]
}
pre <- load_cc("case_crossover_df_timestrat_month_post7.RDS")
cc  <- load_cc("case_crossover_df_shock_post7.RDS")

summ <- function(d) d[, .(cases = sum(outbreak_binary), refs = sum(outbreak_binary == 0),
                          strata = uniqueN(stratum_id), cells = uniqueN(zone_id),
                          per = sum(outbreak_binary == 0) / sum(outbreak_binary)), by = g]
SP <- summ(pre); SC <- summ(cc)
P <- function(grp, f) SP[g == grp][[f]]
C <- function(grp, f) SC[g == grp][[f]]
lost <- function(g) P(g, "cases") - C(g, "cases")

tot_cases  <- sum(cc$outbreak_binary)
tot_refs   <- sum(cc$outbreak_binary == 0)
tot_strata <- uniqueN(cc$stratum_id)
tot_cells  <- uniqueN(cc$zone_id)

# 2a. Layout. Trunk boxes run down the centre; the branch splits at y = 0.50 and rejoins at the
#     bottom. Exclusions hang off to the right of whatever they exclude from.
TX1 <- 0.22; TX2 <- 0.78; TXM <- (TX1 + TX2) / 2           # trunk box left / right / middle
LX1 <- 0.02; LX2 <- 0.47; LXM <- (LX1 + LX2) / 2           # left branch
RX1 <- 0.53; RX2 <- 0.98; RXM <- (RX1 + RX2) / 2           # right branch
HB  <- 0.093

trunk <- data.table(
  x = TXM, xmin = TX1, xmax = TX2, y = c(0.950, 0.805, 0.655, 0.070), h = HB,
  label = c(
    sprintf("Daily 10 km grid, %d northern Mississippi Flyway states\n%s to %s\n%s grid cells, %s cell-days",
            n_states, date_from, date_to, fmt(n_gridcells), fmt(n_celldays)),
    sprintf("USDA APHIS confirmed HPAI poultry detections\n%d events → %d case cell-days in %d cells",
            n_raw_events, n_raw_celldays, n_raw_cells),
    sprintf("Eligible spillover events\n%d events → %d case cell-days in %d cells",
            n_elig_events, n_elig_celldays, n_elig_cells),
    sprintf("Pooled analytic sample\n%d cases  |  %s referent cell-days  |  %d strata in %d cells",
            tot_cases, fmt(tot_refs), tot_strata, tot_cells)))

branch <- data.table(
  x    = c(LXM, RXM, LXM, RXM),
  xmin = c(LX1, RX1, LX1, RX1),
  xmax = c(LX2, RX2, LX2, RX2),
  y    = c(0.455, 0.455, 0.240, 0.240), h = HB,
  label = c(
    sprintf("Minnesota\n%d events → %d case cell-days in %d cells",
            E("mn", "events"), P("mn", "cases"), P("mn", "cells")),
    sprintf("Other northern Mississippi Flyway states\n%d events → %d case cell-days in %d cells",
            E("other", "events"), P("other", "cases"), P("other", "cells")),
    sprintf("Minnesota analytic sample\n%d cases  |  %s referents  |  %d strata in %d cells\n%.1f referents per case",
            C("mn", "cases"), fmt(C("mn", "refs")), C("mn", "strata"), C("mn", "cells"), C("mn", "per")),
    sprintf("Other flyway analytic sample\n%d cases  |  %s referents  |  %d strata in %d cells\n%.1f referents per case",
            C("other", "cases"), fmt(C("other", "refs")), C("other", "strata"), C("other", "cells"),
            C("other", "per"))))

# 2b. Exclusion boxes, hanging right off the trunk / branches
excl <- data.table(
  x = c(TXM + 0.30, LXM + 0.135, RXM + 0.135), y = c(0.733, 0.348, 0.348),
  xmin = c(TXM + 0.06, LXM + 0.03, RXM + 0.03),
  xmax = c(TXM + 0.54, LXM + 0.24, RXM + 0.24),
  h = c(0.068, 0.060, 0.060),
  label = c(
    sprintf("Excluded: farm-to-farm transmission\nflagged by sequence and epidemiological\ninvestigation  (− %d events)", n_flagged),
    sprintf("Excluded: no 90-day baseline\n(− %d cases)", lost("mn")),
    sprintf("Excluded: no 90-day baseline\n(− %d cases)", lost("other"))))

# 2c. Arrows. Down the trunk, out to the split, and back together at the bottom.
SPLIT <- 0.545; JOIN <- 0.145
seg <- rbind(
  data.table(x = TXM, xend = TXM, y = trunk$y[1] - HB/2, yend = trunk$y[2] + HB/2),
  data.table(x = TXM, xend = TXM, y = trunk$y[2] - HB/2, yend = trunk$y[3] + HB/2),
  data.table(x = LXM, xend = LXM, y = SPLIT, yend = 0.455 + HB/2),
  data.table(x = RXM, xend = RXM, y = SPLIT, yend = 0.455 + HB/2),
  data.table(x = LXM, xend = LXM, y = 0.455 - HB/2, yend = 0.240 + HB/2),
  data.table(x = RXM, xend = RXM, y = 0.455 - HB/2, yend = 0.240 + HB/2),
  data.table(x = LXM, xend = LXM, y = 0.240 - HB/2, yend = JOIN),
  data.table(x = RXM, xend = RXM, y = 0.240 - HB/2, yend = JOIN))
# the split and the rejoin are plain connectors, no arrowheads
plain <- rbind(
  data.table(x = TXM, xend = TXM, y = trunk$y[3] - HB/2, yend = SPLIT),
  data.table(x = LXM, xend = RXM, y = SPLIT, yend = SPLIT),
  data.table(x = LXM, xend = RXM, y = JOIN,  yend = JOIN))
join <- data.table(x = TXM, xend = TXM, y = JOIN, yend = trunk$y[4] + HB/2)
# short stubs from the trunk / branches into each exclusion box
stub <- data.table(x = c(TXM, LXM, RXM), xend = c(TXM + 0.06, LXM + 0.03, RXM + 0.03),
                   y = c(0.733, 0.348, 0.348), yend = c(0.733, 0.348, 0.348))

allbox <- rbind(trunk[, .(x, xmin, xmax, y, h, label)],
                branch[, .(x, xmin, xmax, y, h, label)])

# 3a. Draw it. No title or caption by design; those live in the manuscript text.
fp <- ggplot() +
  geom_rect(data = allbox, aes(xmin = xmin, xmax = xmax, ymin = y - h/2, ymax = y + h/2),
            fill = "white", colour = "grey25", linewidth = 0.4) +
  geom_rect(data = excl, aes(xmin = xmin, xmax = xmax, ymin = y - h/2, ymax = y + h/2),
            fill = "grey96", colour = "grey55", linewidth = 0.35, linetype = "dashed") +
  geom_segment(data = rbind(seg, join), aes(x = x, xend = xend, y = y, yend = yend),
               arrow = arrow(length = unit(0.17, "cm"), type = "closed"),
               colour = "grey25", linewidth = 0.4) +
  geom_segment(data = plain, aes(x = x, xend = xend, y = y, yend = yend),
               colour = "grey25", linewidth = 0.4) +
  geom_segment(data = stub, aes(x = x, xend = xend, y = y, yend = yend),
               colour = "grey55", linewidth = 0.35, linetype = "dashed") +
  geom_text(data = allbox, aes(x = x, y = y, label = label),
            family = "Avenir", size = 2.95, lineheight = 1.14) +
  geom_text(data = excl, aes(x = x, y = y, label = label),
            family = "Avenir", size = 2.5, colour = "grey25", lineheight = 1.14) +
  annotate("text", x = TXM, y = 0.578, family = "Avenir", size = 2.55, colour = "grey30",
           lineheight = 1.15,
           label = "Referents: every other day in the same calendar month within the same cell,\nexcluding the 7 days after detection") +
  scale_x_continuous(limits = c(0, 1), expand = c(0, 0)) +
  scale_y_continuous(limits = c(0, 1), expand = c(0, 0)) +
  theme_spark_map()

# 4a. Export. PDF through cairo_pdf, since the default pdf() device can't embed Avenir.
ggsave_spark(file.path(figures_main_folder, "figure1_consort_flowchart.png"),
             fp, width = 10.5, height = 8)
ggsave(file.path(figures_main_folder, "figure1_consort_flowchart.pdf"),
       fp, width = 10.5, height = 8, bg = "white", device = grDevices::cairo_pdf)

cat(sprintf("%d events -> minus %d flagged -> %d eligible -> MN %d / other %d cases\n",
            n_raw_events, n_flagged, n_elig_events, C("mn", "cases"), C("other", "cases")))
cat("wrote figure1_consort_flowchart.png/.pdf\n")

# 5. The case-cell lookup the map (d_02) draws from. It used to be a hand-made file with no
#    producer; it is exactly the eligible case cells, so it is written here from the same panel
#    the flowchart counts. 161 rows: zone_id, state, layer.
el <- setDT(readRDS(paste0(objects_folder, "case_crossover_df_timestrat_month_post7.RDS")))
if (anyDuplicated(names(el))) el <- el[, which(!duplicated(names(el))), with = FALSE]
ml <- unique(el[outbreak_binary == 1, .(zone_id, state)])[, layer := "case"]
fwrite(ml, paste0(objects_folder, "map_layers_flyway.csv"))
cat(sprintf("wrote map_layers_flyway.csv (%d case cells)\n", nrow(ml)))
