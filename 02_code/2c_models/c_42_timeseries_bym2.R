############ MN TIME-SERIES DLNM v2 — BYM2, ZERO-INFLATION, NO BIRDS, SIMPLER CROSSBASIS ############
# Four fixes to c_41, each addressing a specific thing that was wrong with it.
#
# 1. BYM2 instead of iid. c_41's intervals were not credible - wind speed 1.00 (0.98, 1.01) off
#    86 events - because an iid cell effect handles each cell's LEVEL but not the correlation
#    BETWEEN cells. The exposure is a smooth climatological field, so neighbouring cells on the
#    same day carry almost the same anomaly and the effective sample size is a small fraction of
#    the 1.28M nominal rows. BYM2 (Riebler et al. 2016) splits the cell effect into spatially
#    structured and unstructured parts with an interpretable mixing parameter.
#
# 2. Zero inflation. 1,189 of 3,334 cells have no registered feedlot at all and cannot host a
#    poultry outbreak, so their zeros are STRUCTURAL, not small-mean Poisson zeros. That is a
#    real mechanism, not a fudge. Caveat: INLA's zero-inflated families carry a single constant
#    inflation probability, so p cannot be made a function of feedlot count here.
#
# 3. Birds dropped. In this design the bird term is squeezed from both sides - the cell effect
#    takes the between-cell habitat signal, the seasonal term takes the flyway-wide migration
#    wave - and it collapsed from 2.53 in the case-crossover to 1.03 in c_41. It contributes
#    nothing and costs parameters. The eBird abundance surface is itself modelled and spatially
#    smooth, so it would also compete with BYM2 for the same variation.
#
# 4. Simpler crossbasis. 86 events cannot support ns(2) x ns(2) on four exposures: that is 16
#    crossbasis parameters before the spatial field, the seasonal term and the feedlot slope,
#    about 5 events per parameter. The ladder below is fitted in full and reported whatever it
#    says, rather than picking a winner after the fact.

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

# ---- 1. panel + climatology, joined in space (the two grids share no zone_id scheme) ----
d <- setDT(readRDS(paste0(objects_folder, "fishnet_timeseries_daily_modeling.RDS")))
d[, date := as.Date(date)]
setorder(d, zone_id, date)

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
cat(sprintf("panel: %s cell-days, %d cells, %d outbreaks\n",
            format(nrow(d), big.mark = ","), uniqueN(d$zone_id), sum(d$outbreak_count)))

# ---- 2. neighbour graph ----
# The fishnet is an exactly regular 0.09 deg lattice (66 x 87, 58% filled by Minnesota's
# outline), so adjacency comes from integer row/col indices rather than a distance threshold -
# no arbitrary cutoff, no risk of a diagonal sneaking in at one latitude and not another.
# Queen contiguity: the 8 surrounding cells.
gr <- unique(d[, .(zone_id, lat, lon)])
setorder(gr, lat, lon)
gr[, `:=`(r = as.integer(round((lat - min(lat)) / 0.09)) + 1L,
          cc = as.integer(round((lon - min(lon)) / 0.09)) + 1L)]
gr[, idx := .I]
key <- gr[, setNames(idx, paste(r, cc))]
OFF <- expand.grid(dr = -1:1, dc = -1:1); OFF <- OFF[!(OFF$dr == 0 & OFF$dc == 0), ]
nb <- rbindlist(lapply(seq_len(nrow(OFF)), function(k) {
  j <- key[paste(gr$r + OFF$dr[k], gr$cc + OFF$dc[k])]
  data.table(i = gr$idx[!is.na(j)], j = j[!is.na(j)])
}))
cat(sprintf("adjacency: %d links, mean %.2f neighbours, %d isolated cells\n",
            nrow(nb), nrow(nb) / nrow(gr), sum(!gr$idx %in% nb$i)))

# INLA wants a symmetric sparse adjacency; components matter because BYM2 needs a connected
# graph (it scales the structured variance per component).
A <- sparseMatrix(i = nb$i, j = nb$j, x = 1, dims = c(nrow(gr), nrow(gr)))
A <- pmin(A + t(A), 1)
g <- INLA::inla.read.graph(A)
cat(sprintf("graph: %d nodes, %d connected components\n", g$n, g$cc$n))

d <- merge(d, gr[, .(zone_id, cell = idx)], by = "zone_id", all.x = TRUE)
setorder(d, zone_id, date)
stopifnot(!anyNA(d$cell))

# ---- 3. crossbasis ladder ----
shape <- function(s) if (s == "lin") list(fun = "lin") else list(fun = "ns", df = 2)
build_bases <- function(vs, ls) {
  lapply(setNames(AN, AN), function(v) {
    m <- scale(as.matrix(d[, paste0(v, "_w", 0:(NWEEK - 1)), with = FALSE]))
    cb <- dlnm::crossbasis(m, lag = c(0, NWEEK - 1), argvar = shape(vs), arglag = shape(ls))
    colnames(cb) <- paste0("cb_", v, ".", colnames(cb)); cb
  })
}
PC   <- list(prec = list(prior = "pc.prec", param = c(1, 0.01)))
BYM2 <- list(prec = list(prior = "pc.prec", param = c(1, 0.01)),
             phi  = list(prior = "pc", param = c(0.5, 0.5)))

fit_one <- function(bases, family) {
  X <- do.call(cbind, unname(lapply(bases, function(b) as.data.frame(unclass(b)))))
  md <- data.frame(y = d$outbreak_count, X, check.names = FALSE,
                   cell = d$cell, doy = as.integer(d$doy), lfeed = log1p(d$feedlot_count))
  cbn <- unlist(lapply(bases, colnames), use.names = FALSE)
  stopifnot(all(cbn %in% names(md)))
  f <- as.formula(paste("y ~", paste(cbn, collapse = " + "), "+ lfeed +",
                        "f(cell, model = 'bym2', graph = g, scale.model = TRUE, hyper = BYM2) +",
                        "f(doy, model = 'rw2', cyclic = TRUE, hyper = PC)"))
  INLA::inla(f, family = family, data = md,
             control.fixed = list(prec = 4, prec.intercept = 0.01),
             control.compute = list(waic = TRUE, config = TRUE),
             control.predictor = list(compute = FALSE),
             control.inla = list(strategy = "gaussian", int.strategy = "eb"))
}

rr <- function(fit, bases, v, w) {
  cb <- bases[[v]]
  bs <- function(x, a) do.call(dlnm::onebasis, c(list(x = x), a))
  ev <- as.numeric(bs(0.5, attr(cb, "argvar")) - bs(0, attr(cb, "argvar")))
  lv <- colSums(bs(w, attr(cb, "arglag")), na.rm = TRUE)
  wt <- as.vector(t(outer(ev, lv)))
  cv <- dlnm_coef_vcov(fit, colnames(cb))
  lr <- sum(wt * cv$coef); se <- sqrt(as.numeric(t(wt) %*% cv$vcov %*% wt))
  sprintf("%.2f (%.2f, %.2f)", exp(lr), exp(lr - 1.96 * se), exp(lr + 1.96 * se))
}

# stage 1: choose the basis under Poisson, the cheapest family, and report the whole ladder
LADDER <- list(c("lin", "lin"), c("lin", "ns2"), c("ns2", "ns2"))
lad <- rbindlist(lapply(LADDER, function(s) {
  b <- build_bases(s[1], s[2]); np <- length(unlist(lapply(b, colnames)))
  cat(sprintf("  basis %-3s x %-3s (%2d crossbasis params) ... ", s[1], s[2], np)); flush.console()
  f <- try(fit_one(b, "poisson"), silent = TRUE)
  if (inherits(f, "try-error")) { cat("failed\n"); return(NULL) }
  cat(sprintf("WAIC %.1f\n", f$waic$waic))
  data.table(var_shape = s[1], lag_shape = s[2], params = np, WAIC = round(f$waic$waic, 1))
}))
cat("\n=== crossbasis ladder (Poisson + BYM2) ===\n"); print(lad)
setorder(lad, WAIC); best <- lad[1]
cat(sprintf("selected: %s x %s\n\n", best$var_shape, best$lag_shape))

# stage 2: families, at the selected basis
bases <- build_bases(best$var_shape, best$lag_shape)
FAMS <- c(Poisson = "poisson", `Negative binomial` = "nbinomial",
          `Zero-inflated Poisson` = "zeroinflatedpoisson1",
          `Zero-inflated NB` = "zeroinflatednbinomial1")
out <- list(); hyp <- list()
for (nm in names(FAMS)) {
  cat(sprintf("--- %s ---\n", nm)); flush.console()
  f <- try(fit_one(bases, FAMS[[nm]]), silent = TRUE)
  if (inherits(f, "try-error")) { cat("failed:", attr(f, "condition")$message, "\n"); next }
  cat(sprintf("WAIC %.1f\n", f$waic$waic))
  h <- f$summary.hyperpar
  print(round(h[, c("mean", "0.025quant", "0.975quant")], 3))
  hyp[[nm]] <- data.table(family = nm, param = rownames(h),
                          mean = h$mean, lo = h$`0.025quant`, hi = h$`0.975quant`)
  out[[nm]] <- rbindlist(lapply(names(bases), function(v) rbindlist(lapply(0:(NWEEK - 1), function(w)
    data.table(Family = nm, WAIC = round(f$waic$waic, 1),
               Exposure = unname(LAB[sub("_an$", "", v)]),
               `Lag week` = paste0("Week ", w), est = rr(f, bases, v, w))))))
  cat("\n")
}
res <- rbindlist(out)
saveRDS(list(ladder = lad, estimates = res, hyper = rbindlist(hyp)),
        paste0(objects_folder, "si_timeseries_bym2.RDS"))
cat("=== Minnesota time-series DLNM, BYM2, 25-year absolute anomaly ===\n")
print(dcast(res, Exposure + Family ~ `Lag week`, value.var = "est"), nrows = 100)
