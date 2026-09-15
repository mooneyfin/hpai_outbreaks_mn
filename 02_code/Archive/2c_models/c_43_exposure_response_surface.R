############ EXPOSURE-RESPONSE SURFACE FROM THE BYM2 TIME-SERIES FIT ############
# c_42 reported a +0.5 SD LINEAR CONTRAST for every exposure and every one came back ~1.00. But
# the same fit chose ns(2) x ns(2) over the linear bases by 107 WAIC points, and a flat surface
# cannot do that. A non-monotone response - harm at one end, protection at the other - averages
# to nothing at a single symmetric contrast, so "null at +0.5 SD" is all c_42 actually showed.
#
# This pulls the whole surface out instead of one number off it.
#
# Two things the c_42 summary got loose and this fixes:
#  - the reference point. The lag matrix went through scale(), so z = 0 is the MEAN anomaly over
#    the panel, not "conditions equal to the 25-year normal". Those differ. Everything here is
#    centred on raw anomaly = 0, which is the interpretable reference.
#  - the units. Curves are drawn in natural units (degC, mm, m3/m3, m/s), not SD.
#
# The fit and its bases are cached this time so nothing downstream has to refit.

project.folder = paste0(print(here::here()), '/')
src <- readLines(paste0(project.folder, 'create_folder_structure.R'))
eval(parse(text = paste(src[!grepl("^\\s*rm\\(list", src)], collapse = "\n")))
source(paste0(functions.folder, 'script_initiate.R'))
suppressMessages({library(INLA); library(FNN); library(Matrix)})
source(paste0(functions.folder, 'inla_dlnm_helpers.R'))
INLA::inla.setOption(num.threads = 4)

NWEEK <- 4
VARS  <- c("temperature", "precipitation", "soil_moisture", "wind_speed")
LAB   <- setNames(c("Temperature", "Precipitation", "Soil moisture", "Wind speed"), VARS)
# ERA5 ships Kelvin and metres; the anomaly is a DIFFERENCE so Kelvin and Celsius coincide
UNIT  <- setNames(c("°C", "mm", "m³/m³", "m/s"), VARS)
MULT  <- setNames(c(1, 1000, 1, 1), VARS)     # m -> mm for precipitation
FIT_RDS <- paste0(objects_folder, "timeseries_bym2_poisson_fit.RDS")

# ---- 1. panel, identical to c_42 ----
d <- setDT(readRDS(paste0(objects_folder, "fishnet_timeseries_daily_modeling.RDS")))
d[, date := as.Date(date)]; setorder(d, zone_id, date)
MEAN <- paste0(VARS, "_mean")
cfile <- list.files(env_data, pattern = "flyway8_climatology_1997_2021_w", full.names = TRUE)
clim <- rbindlist(lapply(cfile, fread, select = c("zone_id", "state", "week_idx", MEAN),
                         colClasses = list(character = "zone_id")))[state == "Minnesota"][, state := NULL]
hubxy <- unique(rbindlist(lapply(cfile, fread, select = c("zone_id", "state", "lat", "lon"),
                                 colClasses = list(character = "zone_id")))[
                                   state == "Minnesota", .(zone_id, lat, lon)])
grid <- unique(d[, .(zone_id, lat, lon)])
kx <- cos(mean(grid$lat) * pi / 180)
nn <- get.knnx(cbind(hubxy$lon * kx, hubxy$lat), cbind(grid$lon * kx, grid$lat), k = 1)
grid[, clim_id := hubxy$zone_id[nn$nn.index]]
d <- merge(d, grid[, .(zone_id, clim_id)], by = "zone_id", all.x = TRUE)
d[, week_idx := pmin(as.integer(strftime(date, "%j")) %/% 7L, 51L)]
d <- merge(d, clim, by.x = c("clim_id", "week_idx"), by.y = c("zone_id", "week_idx"), all.x = TRUE)
setorder(d, zone_id, date)
for (v in VARS) set(d, j = paste0(v, "_an"), value = d[[v]] - d[[paste0(v, "_mean")]])
AN <- paste0(VARS, "_an")
roll7 <- function(x) data.table::frollmean(x, 7, align = "right")
for (v in AN) {
  d[, (paste0(v, "_w0")) := roll7(get(v)), by = zone_id]
  for (w in 1:(NWEEK - 1)) d[, (paste0(v, "_w", w)) := shift(get(paste0(v, "_w0")), 7L * w), by = zone_id]
}
need <- as.vector(outer(AN, paste0("_w", 0:(NWEEK - 1)), paste0))
d <- d[complete.cases(d[, ..need])]

gr <- unique(d[, .(zone_id, lat, lon)]); setorder(gr, lat, lon)
gr[, `:=`(r = as.integer(round((lat - min(lat)) / 0.09)) + 1L,
          cc = as.integer(round((lon - min(lon)) / 0.09)) + 1L)][, idx := .I]
key <- gr[, setNames(idx, paste(r, cc))]
OFF <- expand.grid(dr = -1:1, dc = -1:1); OFF <- OFF[!(OFF$dr == 0 & OFF$dc == 0), ]
nb <- rbindlist(lapply(seq_len(nrow(OFF)), function(k) {
  j <- key[paste(gr$r + OFF$dr[k], gr$cc + OFF$dc[k])]
  data.table(i = gr$idx[!is.na(j)], j = j[!is.na(j)])
}))
A <- sparseMatrix(i = nb$i, j = nb$j, x = 1, dims = c(nrow(gr), nrow(gr)))
A <- pmin(A + t(A), 1); g <- INLA::inla.read.graph(A)
d <- merge(d, gr[, .(zone_id, cell = idx)], by = "zone_id", all.x = TRUE); setorder(d, zone_id, date)

# ---- 2. bases, keeping the centre/scale so we can get back to natural units ----
bases <- list(); ctr <- numeric(0); scl <- numeric(0)
for (v in AN) {
  m0 <- as.matrix(d[, paste0(v, "_w", 0:(NWEEK - 1)), with = FALSE])
  m  <- scale(m0)
  ctr[v] <- attr(m, "scaled:center")[1]; scl[v] <- attr(m, "scaled:scale")[1]
  cb <- dlnm::crossbasis(m, lag = c(0, NWEEK - 1),
                         argvar = list(fun = "ns", df = 2), arglag = list(fun = "ns", df = 2))
  colnames(cb) <- paste0("cb_", v, ".", colnames(cb))
  bases[[v]] <- cb
}

# ---- 3. fit once, cache ----
if (file.exists(FIT_RDS)) {
  cat("loading cached fit\n"); fit <- readRDS(FIT_RDS)
} else {
  X <- do.call(cbind, unname(lapply(bases, function(b) as.data.frame(unclass(b)))))
  md <- data.frame(y = d$outbreak_count, X, check.names = FALSE,
                   cell = d$cell, doy = as.integer(d$doy), lfeed = log1p(d$feedlot_count))
  cbn <- unlist(lapply(bases, colnames), use.names = FALSE)
  stopifnot(all(cbn %in% names(md)))
  PC   <- list(prec = list(prior = "pc.prec", param = c(1, 0.01)))
  BYM2 <- list(prec = list(prior = "pc.prec", param = c(1, 0.01)),
               phi  = list(prior = "pc", param = c(0.5, 0.5)))
  f <- as.formula(paste("y ~", paste(cbn, collapse = " + "), "+ lfeed +",
                        "f(cell, model = 'bym2', graph = g, scale.model = TRUE, hyper = BYM2) +",
                        "f(doy, model = 'rw2', cyclic = TRUE, hyper = PC)"))
  cat("fitting Poisson + BYM2 ...\n"); flush.console()
  fit <- INLA::inla(f, family = "poisson", data = md,
                    control.fixed = list(prec = 4, prec.intercept = 0.01),
                    control.compute = list(waic = TRUE, config = TRUE),
                    control.predictor = list(compute = FALSE),
                    control.inla = list(strategy = "gaussian", int.strategy = "eb"))
  saveRDS(fit, FIT_RDS)
}
cat(sprintf("WAIC %.1f\n\n", fit$waic$waic))

# ---- 4. predict the surface ----
# grid spans the 1st-99th percentile of each observed lag-0 anomaly, so nothing is extrapolated
pred <- list(); curves <- list(); surf <- list()
for (v in AN) {
  x   <- d[[paste0(v, "_w0")]]
  q   <- quantile(x, c(0.01, 0.99), na.rm = TRUE)
  raw <- seq(q[1], q[2], length.out = 60)
  z   <- (raw - ctr[v]) / scl[v]
  cen_z <- (0 - ctr[v]) / scl[v]                 # reference = equal to the 25-year normal
  cv  <- dlnm_coef_vcov(fit, colnames(bases[[v]]))
  cp  <- dlnm::crosspred(bases[[v]], coef = cv$coef, vcov = cv$vcov,
                         model.link = "log", at = z, cen = cen_z, cumul = TRUE)
  vv  <- sub("_an$", "", v); u <- MULT[[vv]]
  curves[[v]] <- data.table(exposure = unname(LAB[vv]), unit = unname(UNIT[vv]),
                            x = raw * u, RR = cp$allRRfit,
                            low = cp$allRRlow, high = cp$allRRhigh)
  surf[[v]] <- rbindlist(lapply(seq_len(NWEEK) - 1L, function(w) data.table(
    exposure = unname(LAB[vv]), lag_week = w, x = raw * u,
    RR = cp$matRRfit[, w + 1], low = cp$matRRlow[, w + 1], high = cp$matRRhigh[, w + 1])))
}
cur <- rbindlist(curves); sf <- rbindlist(surf)
saveRDS(list(curves = cur, surface = sf, centre = ctr, scale = scl),
        paste0(objects_folder, "si_exposure_response_surface.RDS"))

# ---- 5. does the cumulative curve leave 1 anywhere? ----
cat("=== overall cumulative exposure-response, where the 95% CrI excludes 1 ===\n")
sig <- cur[low > 1 | high < 1]
if (!nrow(sig)) cat("nowhere, for any exposure, across the 1st-99th percentile\n") else
  print(sig[, .(exposure, x = round(x, 3), est = sprintf("%.2f (%.2f, %.2f)", RR, low, high))])
cat("\n=== range of the cumulative RR across the observed exposure range ===\n")
print(cur[, .(unit = unit[1], min_x = round(min(x), 2), max_x = round(max(x), 2),
              RR_min = round(min(RR), 2), RR_max = round(max(RR), 2),
              widest = sprintf("%.2f (%.2f, %.2f)", RR[which.max(abs(log(RR)))],
                               low[which.max(abs(log(RR)))], high[which.max(abs(log(RR)))])),
          by = exposure])

# ---- 6. figure ----
cur[, exposure := factor(exposure, levels = unname(LAB))]
xlabs <- cur[, .(lab = sprintf("%s (%s)", exposure[1], unit[1])), by = exposure]
fig <- ggplot(cur, aes(x, RR)) +
  geom_hline(yintercept = 1, colour = "grey45", linewidth = 0.3) +
  geom_vline(xintercept = 0, colour = "grey75", linewidth = 0.3, linetype = "dashed") +
  geom_ribbon(aes(ymin = low, ymax = high), fill = "#31688E", alpha = 0.16) +
  geom_line(colour = "#31688E", linewidth = 0.8) +
  facet_wrap(~ exposure, nrow = 1, scales = "free_x") +
  scale_y_continuous(trans = "log", breaks = c(0.5, 0.75, 1, 1.5, 2, 3),
                     labels = scales::label_number(drop0trailing = TRUE)) +
  labs(x = "Departure from the 25-year weekly normal (natural units)",
       y = "Cumulative rate ratio over lag weeks 0–3") +
  theme_classic(base_family = "Avenir") +
  theme(strip.text = element_text(family = "Avenir Heavy", size = 11),
        strip.background = element_blank(),
        axis.text = element_text(size = 9.5, colour = "black"),
        axis.title = element_text(size = 11),
        panel.border = element_rect(fill = NA, colour = "grey20", linewidth = 0.4))
ggsave_spark(file.path(figures_main_folder, "figureS3_exposure_response_surface.png"),
             fig, width = 13, height = 4)
cat("\nwrote figureS3_exposure_response_surface.png\n")
