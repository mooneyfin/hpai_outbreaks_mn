############ THE INTERACTION, ONE PRE-STATED SPECIFICATION ############
# c_46 to c_49 explored the specification space. This is the single model the paper reports.
# Every choice below is made on methodological grounds and stated BEFORE the estimates, so none
# of it is picked on the answer it gives.
#
#   exposure     BOTH absolute conditions and the 90-day shock, reported side by side. That is
#                how the primary results are presented everywhere else in the paper, and the two
#                answer different questions: absolute asks whether the conditions themselves
#                mattered, the shock asks whether a CHANGE from recent conditions did. The
#                hypothesis is about "hits", so the shock is the mechanism-matched one, but
#                showing only that would be choosing the parameterisation.
#   lag basis    strata with a break at day 8. The hypothesis names two windows, 0-7 and 8-14,
#                and strata estimates them as independent piecewise-constant bins. A spline
#                basis smooths across the break, so a real effect in one window mechanically
#                drags the other; that is not a shape we want imposed on a windowed hypothesis.
#   modifier     continuous, not dichotomised. A median or tertile split discards information
#                and invents a cut point. The continuous product term is the actual interaction.
#   bird window  lag 15-28 d (weeks 3-4), weather lag 0-14 d (weeks 1-2): the two windows do
#                not overlap, so a day's weather is never also part of the bird term that
#                modifies it. 21-28, 14-28 and 8-28 d are reported as sensitivities (Table S7b).
#   met windows  0-7 d and 8-14 d, pre-stated ("the two weeks prior").
#   prior        INLA_PREC_PRIMARY, the same prior as the paper's primary model. The exploratory
#                scripts used prec = 4, which was an inconsistency, not a decision.
#   panels       all three, as everywhere else in the paper.
#
# REPORTING. Rate ratio with a 95% credible interval AND the posterior probability of direction.
# P(dir) is reported as a continuous quantity, not thresholded: it lets a 0.97 read as "probably
# real, not established" instead of being binned as null. It is computed by sampling the joint
# posterior rather than from a normal approximation to it.

project.folder = paste0(print(here::here()), '/')
src <- readLines(paste0(project.folder, 'create_folder_structure.R'))
eval(parse(text = paste(src[!grepl("^\\s*rm\\(list", src)], collapse = "\n")))
source(paste0(functions.folder, 'script_initiate.R'))
suppressMessages({library(INLA); library(Matrix)})
INLA::inla.setOption(num.threads = 4)
source(paste0(functions.folder, 'inla_dlnm_helpers.R'))

PREC    <- INLA_PREC_PRIMARY          # same prior as the primary model
MAXLAG  <- 14
# the two pre-stated bins, plus their sum. With a strata basis the bins are independent
# coefficients, so 0-14 is just w_a + w_b read off the same posterior - the interval comes from
# the joint draws, which is why it is computed here rather than multiplied by hand
WINDOW  <- list(`0-7 days` = 0:7, `8-14 days` = 8:14, `0-14 days` = 0:14)
BIRDLAG <- 15:28
ARGLAG  <- list(fun = "strata", breaks = 8)
# exposures come from INLA_PRIMARY_* in inla_dlnm_helpers.R - never redefined locally
VARS    <- INLA_PRIMARY_MET
LAB     <- INLA_PRIMARY_LAB
NSAMP   <- 4000

dat <- setDT(readRDS(paste0(objects_folder, "case_crossover_df_shock_post7.RDS")))
if (anyDuplicated(names(dat))) dat <- dat[, which(!duplicated(names(dat))), with = FALSE]
bl <- paste0("anseriformes_Lag", BIRDLAG)
dat[, bird_dep := rowMeans(.SD, na.rm = TRUE), .SDcols = bl]
dat[, bird_c := bird_dep - mean(bird_dep, na.rm = TRUE), by = stratum_id]
dat[, bird_z := bird_c / sd(bird_c, na.rm = TRUE)]

scale_cb <- function(cb, z) { m <- unclass(cb) * z; attributes(m) <- attributes(cb)
                              dimnames(m) <- dimnames(cb); m }
lagw <- function(cb, lags) {
  bs <- function(x, a) do.call(dlnm::onebasis, c(list(x = x), a))
  ev <- as.numeric(bs(0.5, attr(cb, "argvar")) - bs(0, attr(cb, "argvar")))
  as.vector(t(outer(ev, colSums(bs(lags, attr(cb, "arglag")), na.rm = TRUE))))
}

fit_panel <- function(d, panel, kind) {
  mvars <- if (kind == "90-day shock") paste0(VARS, "_shock") else VARS
  bases <- list(); fr <- data.frame(outbreak_binary = d$outbreak_binary)
  add <- function(nm, m) { bases[[nm]] <<- m; fr <<- cbind(fr, as.data.frame(unclass(m))) }
  for (v in mvars) {
    cb <- build_met_crossbasis(d, v, MAXLAG, list(fun = "ns", df = 2), ARGLAG)
    m1 <- cb; colnames(m1) <- paste0("cb_", v, "_main.", colnames(cb)); add(paste0(v, "_main"), m1)
    m2 <- scale_cb(cb, d$bird_z); colnames(m2) <- paste0("cb_", v, "_int.", colnames(cb))
    add(paste0(v, "_int"), m2)
  }
  cbb <- build_met_crossbasis(d, "anseriformes", 28, list(fun = "lin"), list(fun = "ns", df = 2))
  colnames(cbb) <- paste0("cb_anseriformes.", colnames(cbb)); add("anseriformes", cbb)
  fr$bird_z <- d$bird_z
  fr$id_stratum <- as.integer(as.factor(d$stratum_id))
  obj <- list(data = fr, bases = bases, cb_cols = lapply(bases, colnames), met_vars = names(bases),
              contrast_steps = list(), weekly = list(), weekly_vars = character(0),
              max_lag = MAXLAG, max_lag_by_var = list(), argvar = list(fun = "ns", df = 2),
              argvar_by_var = list(), arglag = ARGLAG, arglag_by_var = list())
  fit <- suppressWarnings(fit_cc_inla_dlnm(obj, fixed_prec = PREC))

  # draw the joint posterior once per panel and reuse it for every contrast
  set.seed(20220101)   # seeded so reruns reproduce the reported digits
  samp <- INLA::inla.posterior.sample(NSAMP, fit, selection = list(), verbose = FALSE)
  nm_all <- rownames(samp[[1]]$latent)
  getdraws <- function(cols) {
    idx <- match(paste0(cols, ":1"), nm_all)
    if (anyNA(idx)) idx <- match(cols, sub(":[0-9]+$", "", nm_all))
    t(vapply(samp, function(s) s$latent[idx, 1], numeric(length(cols))))
  }
  summ <- function(draws) {
    q <- quantile(draws, c(.5, .025, .975))
    list(est = sprintf("%.2f (%.2f, %.2f)", exp(q[1]), exp(q[2]), exp(q[3])),
         pdir = max(mean(draws > 0), mean(draws < 0)))
  }
  rbindlist(lapply(mvars, function(v) rbindlist(lapply(names(WINDOW), function(win) {
    lags <- WINDOW[[win]]
    wt   <- lagw(bases[[paste0(v, "_main")]], lags)
    dm   <- getdraws(colnames(bases[[paste0(v, "_main")]])) %*% wt   # effect at mean birds
    di   <- getdraws(colnames(bases[[paste0(v, "_int")]]))  %*% wt   # per +1 SD of birds
    hi   <- summ(dm + di); lo <- summ(dm - di); rt <- summ(2 * di)
    data.table(Panel = panel, Kind = kind, n = sum(d$outbreak_binary),
               Exposure = unname(LAB[sub("_shock$", "", v)]), Window = win,
               `High waterfowl` = hi$est, `P(dir) high` = round(hi$pdir, 3),
               `Low waterfowl`  = lo$est, `P(dir) low`  = round(lo$pdir, 3),
               `Ratio` = rt$est, `P(dir) ratio` = round(rt$pdir, 3))
  }))))
}

res <- rbindlist(lapply(c("Absolute conditions", "90-day shock"), function(kd)
  rbindlist(lapply(c("Minnesota", "Other flyway states", "Pooled"), function(pn) {
    d <- if (pn == "Minnesota") dat[state == "Minnesota"]
         else if (pn == "Pooled") copy(dat) else dat[state != "Minnesota"]
    cat(sprintf("%-20s %-22s n = %d\n", kd, pn, sum(d$outbreak_binary))); flush.console()
    fit_panel(d, pn, kd)
  }))))
saveRDS(res, paste0(objects_folder, "si_interaction_primary.RDS"))

cat("\n=== meteorological shock by prior waterfowl, 90-day shock, strata lag basis ===\n")
cat("    High/Low = effect at +1 / -1 SD of Anseriformes at lag 21-28 d\n")
cat("    Ratio    = high vs low; P(dir) = posterior probability of direction\n\n")
print(res[, .(Kind, Panel, Exposure, Window, `High waterfowl`, `P(dir) high`,
              `Low waterfowl`, Ratio, `P(dir) ratio`)], nrows = 200)
cat("\n=== lag 8-14 days only, the hypothesised window ===\n")
print(res[Window == "8-14 days", .(Kind, Panel, Exposure, `High waterfowl`, `P(dir) high`,
                                   Ratio, `P(dir) ratio`)], nrows = 100)
