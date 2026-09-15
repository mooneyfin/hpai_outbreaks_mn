############ SPACE-TIME INTERACTION — DOES THE EXPOSURE EFFECT SURVIVE IT? ############
# Section 28 said BYM2 did not widen the intervals and blamed the lack of a space-time
# interaction. This tests that claim rather than asserting it.
#
# WHY BYM2 DID NOTHING, AND WHY AN INTERACTION IS A DIFFERENT ANIMAL
# BYM2 is a random effect on the CELL. It says some cells have more outbreaks than others and
# that neighbouring cells are alike. It is constant in time, so it cannot induce correlation
# between two cell-days that differ only in date - which is where the exposure lives. A
# space-time interaction varies by cell AND by time, so for the first time the model has a term
# that could compete with the exposure for the same variation.
#
# WHAT IS ACTUALLY FEASIBLE
# A full Knorr-Held type IV (structured space x structured time) is cell x day = 1,373,608
# random effects. Not happening here. The tractable rungs:
#     cell x week   196,706   too fine - the exposure is a 4-week lag distribution, so a weekly
#                             cell deviation would absorb the exposure outright
#     cell x month   46,676   fitted below
#
# READ IT AS A BOUND, NOT AS THE ANSWER. A cell-month effect operates at roughly the same
# timescale as the 4-week crossbasis, so this is close to over-adjustment by construction. If an
# effect SURVIVES it, that is strong evidence. If it vanishes, the honest conclusion is that this
# design cannot separate the two, not that the exposure does nothing.

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

# ---- 1. panel + graph, same construction as c_42/c_43 ----
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

d[, ym := as.integer(as.factor(format(date, "%Y-%m")))]
d[, cm := as.integer(as.factor(paste(cell, ym)))]
cat(sprintf("panel %s rows | %d cells | %d months | %d cell-months | %d outbreaks\n",
            format(nrow(d), big.mark = ","), uniqueN(d$cell), uniqueN(d$ym),
            uniqueN(d$cm), sum(d$outbreak_count)))

# ---- 2. bases ----
bases <- lapply(setNames(AN, AN), function(v) {
  m <- scale(as.matrix(d[, paste0(v, "_w", 0:(NWEEK - 1)), with = FALSE]))
  cb <- dlnm::crossbasis(m, lag = c(0, NWEEK - 1),
                         argvar = list(fun = "ns", df = 2), arglag = list(fun = "ns", df = 2))
  colnames(cb) <- paste0("cb_", v, ".", colnames(cb)); cb
})
X   <- do.call(cbind, unname(lapply(bases, function(b) as.data.frame(unclass(b)))))
cbn <- unlist(lapply(bases, colnames), use.names = FALSE)
md  <- data.frame(y = d$outbreak_count, X, check.names = FALSE, cell = d$cell,
                  doy = as.integer(d$doy), cm = d$cm, lfeed = log1p(d$feedlot_count))
stopifnot(all(cbn %in% names(md)))

PC   <- list(prec = list(prior = "pc.prec", param = c(1, 0.01)))
BYM2 <- list(prec = list(prior = "pc.prec", param = c(1, 0.01)),
             phi  = list(prior = "pc", param = c(0.5, 0.5)))
BASE <- paste("f(cell, model = 'bym2', graph = g, scale.model = TRUE, hyper = BYM2)",
              "f(doy, model = 'rw2', cyclic = TRUE, hyper = PC)", sep = " + ")

fit_it <- function(extra) {
  f <- as.formula(paste("y ~", paste(cbn, collapse = " + "), "+ lfeed +", BASE,
                        if (nzchar(extra)) paste("+", extra) else ""))
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
  data.table(RR = exp(lr), low = exp(lr - 1.96 * se), high = exp(lr + 1.96 * se), se = se)
}

MODS <- list(`No interaction (c_42)` = "",
             `Cell x month, iid`     = "f(cm, model = 'iid', hyper = PC)")
out <- rbindlist(lapply(names(MODS), function(nm) {
  cat(sprintf("\n--- %s ---\n", nm)); flush.console()
  f <- try(fit_it(MODS[[nm]]), silent = TRUE)
  if (inherits(f, "try-error")) { cat("failed:", attr(f, "condition")$message, "\n"); return(NULL) }
  cat(sprintf("WAIC %.1f\n", f$waic$waic))
  print(round(f$summary.hyperpar[, c("mean", "0.025quant", "0.975quant")], 3))
  rbindlist(lapply(names(bases), function(v) rbindlist(lapply(0:(NWEEK - 1), function(w)
    cbind(data.table(Model = nm, WAIC = round(f$waic$waic, 1),
                     Exposure = unname(LAB[sub("_an$", "", v)]),
                     `Lag week` = paste0("Week ", w)), rr(f, v, w))))))
}))
out[, est := sprintf("%.2f (%.2f, %.2f)", RR, low, high)]
saveRDS(out, paste0(objects_folder, "si_spacetime_interaction.RDS"))

cat("\n=== rate ratios ===\n")
print(dcast(out, Exposure + Model ~ `Lag week`, value.var = "est"), nrows = 60)
cat("\n=== did the intervals widen? (ratio of SEs, interaction / none) ===\n")
w <- dcast(out, Exposure + `Lag week` ~ Model, value.var = "se")
setnames(w, c("No interaction (c_42)", "Cell x month, iid"), c("se_none", "se_int"), skip_absent = TRUE)
if (all(c("se_none", "se_int") %in% names(w)))
  print(w[, .(Exposure, `Lag week`, se_none = round(se_none, 4),
              se_int = round(se_int, 4), ratio = round(se_int / se_none, 2))])
