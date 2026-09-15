############ IS THE BIRD x METEOROLOGY INTERACTION REAL? TWO CHECKS ############
# c_46/c_50 found the meteorological effect only appears where waterfowl were abundant two to four weeks
# earlier. Before anyone believes it, two things have to be ruled out.
#
# CHECK A - THE LAG BASIS MIGHT BE MANUFACTURING IT.
# c_46 used ns(df = 2) across a 0-14 day lag axis. That is a smooth two-parameter curve, so a real
# effect at 8-14 days mechanically drags the 0-7 estimate the other way - which is exactly the sign
# flip we saw. Refit with arglag = strata(break at 8): two piecewise-constant bins, 0-7 and 8-14,
# estimated independently with no smoothing between them. If 8-14 survives that, the shape was not
# imposed. ns(3) is fitted too as a middle rung.
#
# CHECK B - THE MEDIAN SPLIT MIGHT BE DOING THE WORK.
# Dichotomising at the stratum median is interpretable but arbitrary. Refit with (i) tertiles, top
# third against bottom third with the middle dropped, and (ii) no dichotomisation at all - the
# meteorological crossbasis multiplied by the continuous standardised bird value, which is the
# proper interaction and throws nothing away.

project.folder = paste0(print(here::here()), '/')
src <- readLines(paste0(project.folder, 'create_folder_structure.R'))
eval(parse(text = paste(src[!grepl("^\\s*rm\\(list", src)], collapse = "\n")))
source(paste0(functions.folder, 'script_initiate.R'))
suppressMessages({library(INLA); library(Matrix)})
INLA::inla.setOption(num.threads = 4)
source(paste0(functions.folder, 'inla_dlnm_helpers.R'))

PREC <- INLA_PREC_PRIMARY; MAXLAG <- 14   # same prior as the primary
WINDOW <- list(`0-7 days` = 0:7, `8-14 days` = 8:14)
BIRDLAG <- 15:28   # the primary modifier window (c_50)
# exposures come from INLA_PRIMARY_* in inla_dlnm_helpers.R - never redefined locally
VARS <- INLA_PRIMARY_MET
LAB  <- INLA_PRIMARY_LAB

dat <- setDT(readRDS(paste0(objects_folder, "case_crossover_df_shock_post7.RDS")))
if (anyDuplicated(names(dat))) dat <- dat[, which(!duplicated(names(dat))), with = FALSE]
bl <- paste0("anseriformes_Lag", BIRDLAG)
dat[, bird_dep := rowMeans(.SD, na.rm = TRUE), .SDcols = bl]
dat[, bird_c := bird_dep - mean(bird_dep, na.rm = TRUE), by = stratum_id]   # within-stratum
dat[, bird_z := bird_c / sd(bird_c, na.rm = TRUE)]
dat[, birdhigh := as.integer(bird_dep > median(bird_dep, na.rm = TRUE)), by = stratum_id]
# rank-based tertiles, not cut() on quantiles: small strata have tied bird values, which makes
# the quantile breaks non-unique and cut() throws
dat[, brank := data.table::frank(bird_dep, ties.method = "average") / .N, by = stratum_id]
dat[, btert := fifelse(brank <= 1/3, 1L, fifelse(brank > 2/3, 3L, 2L))]

mask <- function(cb, keep) { m <- unclass(cb) * keep; attributes(m) <- attributes(cb)
                             dimnames(m) <- dimnames(cb); m }
scale_cb <- function(cb, z) { m <- unclass(cb) * z; attributes(m) <- attributes(cb)
                              dimnames(m) <- dimnames(cb); m }
lagw <- function(cb, lags) {
  bs <- function(x, a) do.call(dlnm::onebasis, c(list(x = x), a))
  ev <- as.numeric(bs(0.5, attr(cb, "argvar")) - bs(0, attr(cb, "argvar")))
  as.vector(t(outer(ev, colSums(bs(lags, attr(cb, "arglag")), na.rm = TRUE))))
}
ARGLAG <- list(`ns(2) [c_46]`   = list(fun = "ns", df = 2),
               `strata at 8`    = list(fun = "strata", breaks = 8),
               `ns(3)`          = list(fun = "ns", df = 3))

# assemble whichever design is asked for and read off the contrast
run <- function(d, mvars, arglag, mode, tag, panel, kind) {
  bases <- list(); fr <- data.frame(outbreak_binary = d$outbreak_binary)
  add <- function(nm, m) { bases[[nm]] <<- m; fr <<- cbind(fr, as.data.frame(unclass(m))) }
  for (v in mvars) {
    cb <- build_met_crossbasis(d, v, MAXLAG, list(fun = "ns", df = 2), arglag)
    if (mode == "continuous") {
      m1 <- cb; colnames(m1) <- paste0("cb_", v, "_main.", colnames(cb)); add(paste0(v, "_main"), m1)
      m2 <- scale_cb(cb, d$bird_z); colnames(m2) <- paste0("cb_", v, "_int.", colnames(cb))
      add(paste0(v, "_int"), m2)
    } else {
      for (grp in c("high", "low")) {
        keep <- if (grp == "high") d$grp_hi else 1 - d$grp_hi
        m <- mask(cb, keep); colnames(m) <- paste0("cb_", v, "_", grp, ".", colnames(cb))
        add(paste0(v, "_", grp), m)
      }
    }
  }
  cbb <- build_met_crossbasis(d, "anseriformes", 28, list(fun = "lin"), list(fun = "ns", df = 2))
  colnames(cbb) <- paste0("cb_anseriformes.", colnames(cbb)); add("anseriformes", cbb)
  if (mode == "continuous") fr$bird_z <- d$bird_z
  fr$id_stratum <- as.integer(as.factor(d$stratum_id))
  obj <- list(data = fr, bases = bases, cb_cols = lapply(bases, colnames), met_vars = names(bases),
              contrast_steps = list(), weekly = list(), weekly_vars = character(0),
              max_lag = MAXLAG, max_lag_by_var = list(), argvar = list(fun = "ns", df = 2),
              argvar_by_var = list(), arglag = arglag, arglag_by_var = list())
  fit <- suppressWarnings(fit_cc_inla_dlnm(obj, fixed_prec = PREC))

  rbindlist(lapply(mvars, function(v) rbindlist(lapply(names(WINDOW), function(win) {
    lags <- WINDOW[[win]]
    if (mode == "continuous") {
      cb <- bases[[paste0(v, "_int")]]; wt <- lagw(cb, lags)
      cv <- dlnm_coef_vcov(fit, colnames(cb))
      # ratio of the effect at +1 SD birds to the effect at -1 SD = exp(2 * w'b_int)
      lr <- 2 * sum(wt * cv$coef)
      se <- 2 * sqrt(as.numeric(t(wt) %*% cv$vcov %*% wt))
      hi <- lo <- ""
    } else {
      a <- met_lag_effect(fit, obj, paste0(v, "_high"), lags)
      b <- met_lag_effect(fit, obj, paste0(v, "_low"),  lags)
      hi <- sprintf("%.2f (%.2f, %.2f)", a[["OR"]], a[["low"]], a[["high"]])
      lo <- sprintf("%.2f (%.2f, %.2f)", b[["OR"]], b[["low"]], b[["high"]])
      cols <- c(obj$cb_cols[[paste0(v, "_high")]], obj$cb_cols[[paste0(v, "_low")]])
      cv <- dlnm_coef_vcov(fit, cols); w1 <- lagw(bases[[paste0(v, "_high")]], lags)
      ct <- c(w1, -w1); lr <- sum(ct * cv$coef)
      se <- sqrt(as.numeric(t(ct) %*% cv$vcov %*% ct))
    }
    data.table(Check = tag, Panel = panel, Kind = kind,
               Exposure = unname(LAB[sub("_shock$", "", v)]), Window = win,
               high = hi, low = lo,
               ratio = sprintf("%.2f (%.2f, %.2f)", exp(lr), exp(lr - 1.96 * se), exp(lr + 1.96 * se)),
               p_int = round(2 * pnorm(-abs(lr / se)), 4))
  }))))
}

grid <- CJ(kind = c("Absolute conditions", "90-day shock"),
           panel = c("Minnesota", "Pooled"), sorted = FALSE)
res <- rbindlist(apply(grid, 1, function(g) {
  kind <- g[["kind"]]; panel <- g[["panel"]]
  mvars <- if (kind == "90-day shock") paste0(VARS, "_shock") else VARS
  d0 <- if (panel == "Minnesota") dat[state == "Minnesota"] else copy(dat)
  out <- list()
  # CHECK A: three lag bases, median split throughout
  for (nm in names(ARGLAG)) {
    d <- copy(d0)[, grp_hi := birdhigh]
    cat(sprintf("A %-20s %-10s %-14s\n", kind, panel, nm)); flush.console()
    out[[paste0("A", nm)]] <- run(d, mvars, ARGLAG[[nm]], "split",
                                  paste("Lag basis:", nm), panel, kind)
  }
  # CHECK B: three ways of using the bird variable, ns(2) lag basis throughout
  for (nm in c("Median split", "Tertiles (1 vs 3)", "Continuous (no split)")) {
    d <- copy(d0)
    if (nm == "Median split")      d[, grp_hi := birdhigh]
    if (nm == "Tertiles (1 vs 3)") { d <- d[btert != 2L]; d[, grp_hi := as.integer(btert == 3L)]
                                     d[, `:=`(k = sum(outbreak_binary), n = .N), by = stratum_id]
                                     d <- d[k == 1L & n > 1L] }
    mode <- if (nm == "Continuous (no split)") "continuous" else "split"
    cat(sprintf("B %-20s %-10s %-22s n=%d\n", kind, panel, nm, sum(d$outbreak_binary))); flush.console()
    r <- try(run(d, mvars, list(fun = "ns", df = 2), mode, paste("Split:", nm), panel, kind),
             silent = TRUE)
    if (!inherits(r, "try-error")) out[[paste0("B", nm)]] <- r
  }
  rbindlist(out)
}))
saveRDS(res, paste0(objects_folder, "si_interaction_robustness.RDS"))

cat("\n================ THE HEADLINE: temperature, lag 8-14 days ================\n")
print(res[Exposure == "Temperature" & Window == "8-14 days",
          .(Check, Panel, Kind, high, low, ratio, p_int)], nrows = 60)
cat("\n================ precipitation, lag 8-14 days ================\n")
print(res[Exposure == "Precipitation" & Window == "8-14 days",
          .(Check, Panel, Kind, ratio, p_int)], nrows = 60)
cat(sprintf("\nsignificant interactions overall: %d of %d (expect %.1f by chance)\n",
            sum(res$p_int < 0.05), nrow(res), 0.05 * nrow(res)))
