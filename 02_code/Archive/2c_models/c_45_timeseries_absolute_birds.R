############ TIME-SERIES ON ABSOLUTE CONDITIONS, WITH BIRDS BACK IN ############
# c_42/c_43 ran the 25-year anomaly with birds dropped. This is the other corner: absolute
# conditions, Anseriformes and predators included.
#
# ONE THING TO KEEP IN MIND READING IT
# In a time-series with a cyclic day-of-year term, absolute conditions are NOT the same animal
# they were in the case-crossover. The doy term removes the average seasonal cycle, so what is
# left of absolute temperature is already close to an anomaly. Absolute vs anomaly is mostly a
# distinction without a difference here - which is itself worth showing.
#
# So model B drops the seasonal control entirely. That is a deliberately confounded model, not a
# candidate: temperature and waterfowl both track season, and so does the epidemic. It is here to
# show how much of any "effect" on absolute exposures is just the calendar.

project.folder = paste0(print(here::here()), '/')
src <- readLines(paste0(project.folder, 'create_folder_structure.R'))
eval(parse(text = paste(src[!grepl("^\\s*rm\\(list", src)], collapse = "\n")))
source(paste0(functions.folder, 'script_initiate.R'))
suppressMessages({library(INLA); library(Matrix)})
source(paste0(functions.folder, 'inla_dlnm_helpers.R'))
INLA::inla.setOption(num.threads = 2)

NWEEK <- 4
MET  <- c("temperature", "precipitation", "soil_moisture", "wind_speed")
BIRD <- c("anseriformes", "predators")
EXPO <- c(MET, BIRD)
LAB  <- setNames(c("Temperature", "Precipitation", "Soil moisture", "Wind speed",
                   "Anseriformes", "Predatory birds"), EXPO)

d <- setDT(readRDS(paste0(objects_folder, "fishnet_timeseries_daily_modeling.RDS")))
d[, date := as.Date(date)]; setorder(d, zone_id, date)

roll7 <- function(x) data.table::frollmean(x, 7, align = "right")
for (v in EXPO) {
  d[, (paste0(v, "_w0")) := roll7(get(v)), by = zone_id]
  for (w in 1:(NWEEK - 1)) d[, (paste0(v, "_w", w)) := shift(get(paste0(v, "_w0")), 7L * w), by = zone_id]
}
need <- as.vector(outer(EXPO, paste0("_w", 0:(NWEEK - 1)), paste0))
d <- d[complete.cases(d[, ..need])]
cat(sprintf("panel %s rows | %d cells | %d outbreaks\n",
            format(nrow(d), big.mark = ","), uniqueN(d$zone_id), sum(d$outbreak_count)))

# same queen-contiguity graph as c_42
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

# birds linear, meteorology ns(2) - same shapes the case-crossover uses
bases <- lapply(setNames(EXPO, EXPO), function(v) {
  m <- scale(as.matrix(d[, paste0(v, "_w", 0:(NWEEK - 1)), with = FALSE]))
  cb <- dlnm::crossbasis(m, lag = c(0, NWEEK - 1),
                         argvar = if (v %in% BIRD) list(fun = "lin") else list(fun = "ns", df = 2),
                         arglag = list(fun = "ns", df = 2))
  colnames(cb) <- paste0("cb_", v, ".", colnames(cb)); cb
})
X   <- do.call(cbind, unname(lapply(bases, function(b) as.data.frame(unclass(b)))))
cbn <- unlist(lapply(bases, colnames), use.names = FALSE)
md  <- data.frame(y = d$outbreak_count, X, check.names = FALSE, cell = d$cell,
                  doy = as.integer(d$doy), lfeed = log1p(d$feedlot_count))
stopifnot(all(cbn %in% names(md)))

PC   <- list(prec = list(prior = "pc.prec", param = c(1, 0.01)))
BYM2 <- list(prec = list(prior = "pc.prec", param = c(1, 0.01)),
             phi  = list(prior = "pc", param = c(0.5, 0.5)))
fit_it <- function(seasonal) {
  f <- as.formula(paste("y ~", paste(cbn, collapse = " + "), "+ lfeed +",
    "f(cell, model = 'bym2', graph = g, scale.model = TRUE, hyper = BYM2)",
    if (seasonal) "+ f(doy, model = 'rw2', cyclic = TRUE, hyper = PC)" else ""))
  INLA::inla(f, family = "poisson", data = md,
             control.fixed = list(prec = 4, prec.intercept = 0.01),
             control.compute = list(waic = TRUE, config = TRUE),
             control.predictor = list(compute = FALSE),
             control.inla = list(strategy = "gaussian", int.strategy = "eb"))
}
rr <- function(fit, v, w) {
  cb <- bases[[v]]
  bs <- function(x, a) do.call(dlnm::onebasis, c(list(x = x), a))
  wt <- as.vector(t(outer(as.numeric(bs(0.5, attr(cb, "argvar")) - bs(0, attr(cb, "argvar"))),
                          colSums(bs(w, attr(cb, "arglag")), na.rm = TRUE))))
  cv <- dlnm_coef_vcov(fit, colnames(cb))
  lr <- sum(wt * cv$coef); se <- sqrt(as.numeric(t(wt) %*% cv$vcov %*% wt))
  sprintf("%.2f (%.2f, %.2f)", exp(lr), exp(lr - 1.96 * se), exp(lr + 1.96 * se))
}

MODS <- c(`A. With seasonal control` = TRUE, `B. NO seasonal control (confounded)` = FALSE)
out <- rbindlist(lapply(names(MODS), function(nm) {
  cat(sprintf("\n--- %s ---\n", nm)); flush.console()
  f <- try(fit_it(MODS[[nm]]), silent = TRUE)
  if (inherits(f, "try-error")) { cat("failed\n"); return(NULL) }
  cat(sprintf("WAIC %.1f\n", f$waic$waic))
  rbindlist(lapply(EXPO, function(v) rbindlist(lapply(0:(NWEEK - 1), function(w)
    data.table(Model = nm, WAIC = round(f$waic$waic, 1), Exposure = unname(LAB[v]),
               `Lag week` = paste0("Week ", w), est = rr(f, v, w))))))
}))
saveRDS(out, paste0(objects_folder, "si_timeseries_absolute_birds.RDS"))
cat("\n=== Minnesota time-series, ABSOLUTE conditions + birds ===\n")
print(dcast(out, Exposure + Model ~ `Lag week`, value.var = "est"), nrows = 60)
