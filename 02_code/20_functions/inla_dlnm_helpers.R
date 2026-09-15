############ INLA CASE-CROSSOVER DLNM HELPERS ############
# Bayesian case-crossover DLNM fit in INLA. Conditional-Poisson trick: stratified
# Poisson with a free intercept per case-crossover stratum = conditional-logistic.
# Met exposures enter as fixed-effect crossbasis columns; weekly bird abundance is a
# random-walk-smoothed z-model over weekly lag knots. Dataset-agnostic (MN + multi-state).
# Requires INLA, dlnm, data.table, Matrix.

#1. Exposure sets
# THE canonical exposure set. Every model that reports a primary or sensitivity estimate takes
# its variables from here, so the set cannot drift between scripts. It had drifted badly: c_21,
# c_23 and c_24 carried runoff but no precipitation, while c_33, c_36, c_37 and c_47 carried
# precipitation but no runoff, and nothing flagged it. Change the exposure set HERE and nowhere
# else.
INLA_PRIMARY_MET   <- c("temperature", "precipitation", "soil_moisture", "wind_speed", "runoff")
INLA_PRIMARY_SHOCK <- paste0(INLA_PRIMARY_MET, "_shock")
INLA_PRIMARY_LAB   <- c(temperature = "Temperature", precipitation = "Precipitation",
                        soil_moisture = "Soil moisture", wind_speed = "Wind speed",
                        runoff = "Runoff", anseriformes = "Anseriformes")
# labels keyed by the shock column names too, so sub("_shock$") is never needed at a call site
INLA_PRIMARY_LAB   <- c(INLA_PRIMARY_LAB,
                        setNames(unname(INLA_PRIMARY_LAB[INLA_PRIMARY_MET]), INLA_PRIMARY_SHOCK))

# Panel display names. The model objects key panels as Minnesota / Other flyway / Pooled; the
# paper says "All northern Mississippi Flyway states" for the third, because it is ONE model fit
# on all 175 strata, not a pooling of estimates. Map on the way out, never rename the keys.
INLA_PANEL_LAB <- c(Minnesota = "Minnesota",
                    `Other flyway` = "Other northern Mississippi Flyway states",
                    `Other flyway states` = "Other northern Mississippi Flyway states",
                    Pooled = "All northern Mississippi Flyway states")
norm_panel <- function(x) { y <- unname(INLA_PANEL_LAB[x]); ifelse(is.na(y), x, y) }
INLA_PANEL_ORD <- unname(INLA_PANEL_LAB[c("Minnesota", "Other flyway", "Pooled")])

# the wider pool a_0x builds; not the modelling set
INLA_MET_VARS    <- c("runoff", "soil_moisture", "temperature",
                      "precipitation", "wind_speed", "snow_cover")
# Only anseriformes carries a weekly term in the primary. Predators is dropped: WAIC prefers
# it out (666.1 vs 668.5), it barely moves the waterfowl estimate (1.17 -> 1.19), and
# scavengers feeding on infected carcasses are plausibly a mediator of the waterfowl->farm
# pathway, so adjusting for it over-adjusts the question we care about. It's reported as a
# secondary analysis instead.
INLA_WEEKLY_VARS <- c("anseriformes")
INLA_BIRD_VAR    <- "anseriformes"
INLA_BIRD_KNOTS  <- c(0, 7, 14, 21, 28)              # weekly lag knots (days)

# RW1, not RW2. RW2 only penalises curvature so it happily extrapolates a straight line off
# the end of the lag window, which is where the old "waterfowl peak at week 4" came from;
# RW1 penalises slope and kills that artifact. WAIC agrees (662.4 vs 668.5). A strata
# crossbasis was also tried and is worse on both counts (673.7, much wider intervals).
INLA_BIRD_RW     <- "rw1"
# PC prior on the weekly RW sd: P(sd > 5) = 0.01. The INLA-ish default of P(sd > 1) = 0.01 is
# tight enough to shrink the weekly lag curve to a dead-flat constant, which hides whatever
# lag structure is there. At U=5 a mild gradient appears and it is stable up to U=20, so the
# curve has converged rather than still drifting. WAIC is a wash (1442.4 vs 1441.4).
INLA_BIRD_PC     <- c(5, 0.01)

# Temperature threshold on the natural scale. ERA5-Land ships temperature in KELVIN, and
# each lag column was z-scaled on its own mean/sd, so convert per panel rather than
# hard-coding a z. 10 degC beat 0/5/15/20 on WAIC and keeps ~28% of the lag matrix above it.
INLA_TEMP_THRESHOLD_C <- 10
INLA_TEMP_STEP_C      <- 5    # report temperature per +5 degC above the threshold
temp_threshold_z <- function(df, thr_C = INLA_TEMP_THRESHOLD_C) {
  ((thr_C + 273.15) - mean(df$temperature, na.rm = TRUE)) / sd(df$temperature, na.rm = TRUE)
}

# ns(df=2) lag basis, WAIC-selected. The N(0, 0.25^2) prior on the DLNM coefficients
# (prec 16) is NOT WAIC-selected, whatever this comment used to claim: WAIC actually prefers
# tighter, monotonically, with no interior optimum (905.2 at 0.25^2 vs 902.4 at 0.125^2 on
# the original panel, and the same ordering on every panel since). It's a deliberate
# regularisation choice at EPV ~ 4, where vague priors return implausible ORs. Report it as
# an assumption and show the sensitivity, because the effect sizes do move with it. Everything except temperature is linear on the exposure scale:
# precipitation was tried as threshold and spline and every variant widened its interval
# without changing direction, and it bottoms out at -0.78 so low cuts are degenerate.
INLA_ARGLAG_PRIMARY <- list(fun = "ns", df = 2)
INLA_PREC_PRIMARY   <- 16
INLA_ARGVAR_BY_VAR  <- list()   # empty => assemble_cc_dlnm derives the temperature threshold

#2. VIF among exposures (lag-0, z-scaled). Result: one VIF per predictor.
compute_exposure_vif <- function(df, predictors = c(INLA_MET_VARS, INLA_BIRD_VAR),
                                 lag = 0) {
  df   <- data.table::as.data.table(df)
  cols <- paste0(predictors, "_Lag", lag)
  X    <- as.matrix(df[, ..cols]); colnames(X) <- predictors
  R    <- stats::cor(X, use = "complete.obs")
  vif  <- diag(solve(R))
  list(
    cor = round(R, 2),
    vif = data.table::data.table(Predictor = predictors, VIF = round(vif, 2))[order(-VIF)]
  )
}

#3. Build one crossbasis from the precomputed x_Lag0..x_LagL columns.
build_met_crossbasis <- function(df, var, max_lag = 28,
                                 argvar = list(fun = "lin"),
                                 arglag = INLA_ARGLAG_PRIMARY,
                                 min_nonzero = 0.05) {
  df <- data.table::as.data.table(df)
  X  <- as.matrix(df[, paste0(var, "_Lag", 0:max_lag), with = FALSE])

  # a threshold cut sitting outside the data zeroes the whole basis, which quietly deletes
  # the variable and looks like a great WAIC. bail loudly instead of shipping that.
  if (identical(argvar$fun, "thr")) {
    nz <- if (startsWith(argvar$side, "l")) mean(X < argvar$thr.value, na.rm = TRUE)
          else                              mean(X > argvar$thr.value, na.rm = TRUE)
    if (nz < min_nonzero)
      stop(sprintf("degenerate threshold for %s: cut %.2f (side %s) leaves only %.1f%% of the lag matrix non-zero",
                   var, argvar$thr.value, argvar$side, 100 * nz))
  }
  dlnm::crossbasis(X, lag = c(0, max_lag), argvar = argvar, arglag = arglag)
}

#4. Random-walk structure matrices for the weekly bird z-model (+ ridge for full rank).
rw2_Cmatrix <- function(K, ridge = 1e-5) {   # penalises curvature
  Matrix::Matrix(crossprod(diff(diag(K), differences = 2)) + ridge * diag(K), sparse = TRUE)
}
rw1_Cmatrix <- function(K, ridge = 1e-5) {   # penalises slope (flatter tails)
  Matrix::Matrix(crossprod(diff(diag(K), differences = 1)) + ridge * diag(K), sparse = TRUE)
}

#5. Assemble the model frame: crossbases, weekly RW design matrices, stratum index.
assemble_cc_dlnm <- function(df,
                             met_vars      = INLA_MET_VARS,
                             weekly_vars   = INLA_WEEKLY_VARS,   # each gets its own RW lag smooth
                             bird_knots    = INLA_BIRD_KNOTS,
                             max_lag       = 28,
                             max_lag_by_var = list(),              # per-var lag window override
                             argvar        = list(fun = "lin"),
                             argvar_by_var = INLA_ARGVAR_BY_VAR,   # per-var exposure-response override
                             arglag        = INLA_ARGLAG_PRIMARY,
                             arglag_by_var = list(),               # per-var lag basis override
                             bird_rw       = INLA_BIRD_RW) {
  df      <- data.table::as.data.table(df)
  n       <- nrow(df)
  bird_rw <- match.arg(bird_rw, c("rw1", "rw2"))

  # temperature gets a threshold at INLA_TEMP_THRESHOLD_C, converted to this panel's own
  # z-scale. caller can override by passing temperature in argvar_by_var.
  if (is.null(argvar_by_var[["temperature"]]) && "temperature" %in% met_vars)
    argvar_by_var[["temperature"]] <- list(fun = "thr", side = "h",
                                           thr.value = temp_threshold_z(df))

  # Met crossbases -> named fixed-effect columns cb_<var>.v1.l1 ...
  bases   <- list(); cb_cols <- list()
  dat     <- data.frame(outbreak_binary = df$outbreak_binary)
  for (v in met_vars) {
    av <- argvar_by_var[[v]] %||% argvar
    al <- arglag_by_var[[v]] %||% arglag
    # each exposure can carry its own lag window. a variable whose mechanism runs for days
    # shouldn't be forced onto the same window as one that runs for weeks, and truncating a
    # real effect doesn't just lose it, it leaks into whatever else is still in the model.
    ml <- max_lag_by_var[[v]] %||% max_lag
    cb <- build_met_crossbasis(df, v, max_lag = ml, argvar = av, arglag = al)
    colnames(cb) <- paste0("cb_", v, ".", colnames(cb))
    bases[[v]]   <- cb
    cb_cols[[v]] <- colnames(cb)
    dat <- cbind(dat, as.data.frame(unclass(cb)))
  }

  # One z-model per weekly exposure: exposure at the lag knots, RW-smoothed. Each gets its
  # own index column id_wk{i} so INLA fits them as separate effects.
  Qw     <- if (bird_rw == "rw1") rw1_Cmatrix(length(bird_knots)) else rw2_Cmatrix(length(bird_knots))
  weekly <- list()
  for (i in seq_along(weekly_vars)) {
    v   <- weekly_vars[i]
    Z   <- as.matrix(df[, paste0(v, "_Lag", bird_knots), with = FALSE])
    idx <- paste0("id_wk", i)
    dat[[idx]] <- seq_len(n)
    weekly[[v]] <- list(var = v, Z = Z, Q = Qw, idx = idx, knots = bird_knots)
  }
  dat$id_stratum <- as.integer(as.factor(df$stratum_id))

  # temperature is contrasted per +5 degC, so convert that to this panel's z units;
  # everything else keeps the shared +0.5 SD step
  contrast_steps <- list()
  if ("temperature" %in% met_vars)
    contrast_steps$temperature <- INLA_TEMP_STEP_C / sd(df$temperature, na.rm = TRUE)

  list(
    data = dat, bases = bases, cb_cols = cb_cols, met_vars = met_vars,
    contrast_steps = contrast_steps,
    weekly = weekly, weekly_vars = weekly_vars, bird_var = weekly_vars[1],
    bird_knots = bird_knots, bird_rw = bird_rw,
    max_lag = max_lag, max_lag_by_var = max_lag_by_var,
    argvar = argvar, argvar_by_var = argvar_by_var,
    arglag = arglag, arglag_by_var = arglag_by_var
  )
}

#6. Fit the conditional-Poisson DLNM. fixed_prec = prior precision on DLNM coefficients
#   (16 = N(0, 0.25^2), regularised primary; 0.01 = N(0, 100), vague).
fit_cc_inla_dlnm <- function(obj, verbose = FALSE, fixed_prec = INLA_PREC_PRIMARY,
                             int_strategy = "eb", bird_pc = INLA_BIRD_PC) {
  hyper.iid    <- list(theta = list(prior = "pc.prec", param = bird_pc))
  hyper.strata <- list(prec = list(initial = log(1e-6), fixed = TRUE))   # free stratum intercepts

  cb_names <- unlist(obj$cb_cols)
  # One z-term per weekly exposure; assign its Z / Cmatrix into this env for the formula.
  z_terms <- character(0)
  for (i in seq_along(obj$weekly)) {
    w <- obj$weekly[[i]]
    assign(paste0("wkZ", i), w$Z); assign(paste0("wkQ", i), w$Q)
    z_terms <- c(z_terms, sprintf(
      "f(%s, model = 'z', Z = wkZ%d, Cmatrix = wkQ%d, hyper = hyper.iid)", w$idx, i, i))
  }
  # assemble from a term vector so a model with no met terms (bird-only) still produces a
  # valid formula instead of "... + + f(...)"
  rhs <- c("f(id_stratum, model = 'iid', constr = FALSE, hyper = hyper.strata)",
           cb_names, z_terms)
  form <- paste0("outbreak_binary ~ 1 + ", paste(rhs, collapse = " + "))
  fixed.mean <- as.list(setNames(rep(0,          length(cb_names)), cb_names))
  fixed.prec <- as.list(setNames(rep(fixed_prec, length(cb_names)), cb_names))

  run <- function(ints) INLA::inla(
    formula = stats::as.formula(form), data = obj$data, family = "poisson", verbose = verbose,
    control.fixed = c(list(correlation.matrix = TRUE),   # covariance stored for crosspred
                      if (length(cb_names)) list(mean = fixed.mean, prec = fixed.prec),
                      list(prec.intercept = 0.001, mean.intercept = 0)),
    control.compute = list(dic = TRUE, waic = TRUE, cpo = TRUE, config = TRUE),
    control.inla    = list(strategy = "simplified.laplace", int.strategy = ints,
                           control.vb = list(enable = FALSE)),   # off: harmless divergence from diffuse strata
    control.mode = list(restart = TRUE), num.threads = 1
  )
  fit <- run(int_strategy)

  # On some panels eb/ccd walk into a junk mode where a weekly-smooth precision runs off to
  # ~1e8. That pins the RW flat, so the exposure reads as exactly 1.00 (0.99, 1.01) instead of
  # its real effect, and it usually segfaults the core on the way. Full grid integration finds
  # the right mode. Cheap to check, and it fails loudly rather than shipping a fake null.
  collapsed <- function(f) {
    if (!length(obj$weekly)) return(FALSE)
    p <- f$summary.hyperpar[grep("^Precision for id_wk",
                                 rownames(f$summary.hyperpar)), "mean"]
    length(p) && any(p > 1e6, na.rm = TRUE)
  }
  if (int_strategy != "grid" && collapsed(fit)) {
    warning("weekly smooth collapsed under int.strategy = '", int_strategy,
            "'; refitting on the full grid", call. = FALSE)
    fit <- run("grid")
    # if grid agrees the smooth is flat, that's a real identifiability limit of the panel and
    # not a bad mode, so keep it and say so. The time-stratified design does this genuinely:
    # its within-stratum lag matrix is nearly all level and no tilt (ends uncorrelated at
    # -0.02, 64% of variance in one all-knots-together direction), so there's no shape to find.
    if (collapsed(fit))
      warning("weekly smooth is flat under grid integration too; treat the lag shape as ",
              "unidentified on this panel rather than as a null result", call. = FALSE)
  }

  # separate failure mode, and grid integration does NOT fix it: the smooth can sit at a
  # perfectly ordinary precision and still come back dead flat, which reads as a clean level
  # estimate. fine as a level, but it is NOT a lag-specific result, so say so loudly rather
  # than let a "lag 22-28 d" number get quoted off a flat line.
  if (length(obj$weekly))
    for (v in obj$weekly_vars)
      if (!weekly_shape_identified(fit, obj, v))
        warning(sprintf(paste0("weekly lag smooth for '%s' is flat (knot range below tol): ",
                               "this fit identifies a LEVEL, not a lag shape. do not report ",
                               "lag- or window-specific effects for it."), v), call. = FALSE)
  fit
}

#7. Coef + covariance for one variable's crossbasis (in crossbasis column order).
dlnm_coef_vcov <- function(fit, cb_cols) {
  idx  <- match(cb_cols, fit$names.fixed)
  coef <- fit$summary.fixed[idx, "mean"]
  V    <- fit$misc$lincomb.derived.covariance.matrix
  if (is.null(V)) V <- config_vcov(fit, cb_cols) else V <- V[idx, idx]
  list(coef = coef, vcov = as.matrix(V))
}

#8a. Pick the contrast for a variable. Linear/spline terms get the usual +0.5 SD from the
#    mean. A threshold basis is flat on one side of its cut-point, so at=0.5/cen=0 would
#    sit in the dead zone and hand back OR = 1; centre those at the threshold and step
#    0.5 SD past it instead, which reads as "per +0.5 SD beyond the threshold".
contrast_for <- function(obj, var, step = NULL) {
  step <- step %||% obj$contrast_steps[[var]] %||% 0.5
  av   <- obj$argvar_by_var[[var]]
  if (is.null(av) || !identical(av$fun, "thr")) return(list(at = step, cen = 0))
  thr <- av$thr.value
  if (startsWith(av$side, "l")) list(at = thr - step, cen = thr)
  else                          list(at = thr + step, cen = thr)
}

#8. crosspred for one variable. at/cen default to contrast_for(); pass a grid to `at` to
#   trace a whole exposure-response curve.
crosspred_cc <- function(fit, obj, var, at = NULL, cen = NULL, cumul = TRUE, ...) {
  ct <- contrast_for(obj, var)
  if (is.null(at))  at  <- ct$at
  if (is.null(cen)) cen <- ct$cen
  cv <- dlnm_coef_vcov(fit, obj$cb_cols[[var]])
  dlnm::crosspred(obj$bases[[var]], coef = cv$coef, vcov = cv$vcov,
                  model.link = "log", at = at, cen = cen, cumul = cumul, ...)
}

#9. Weekly RW lag-response for one weekly exposure: OR per +0.5 SD at each knot + cumulative.
#   `window` = c(from, to) days restricts the cumulative to those knots (e.g. c(21,28) for the
#   ~4-week-prior effect); NULL uses all knots (full 0-28 d cumulative).
weekly_lag_summary <- function(fit, obj, var = obj$weekly_vars[1], at = 0.5, n_sample = 2000,
                               window = NULL) {
  w    <- obj$weekly[[var]]; K <- length(w$knots)
  b    <- utils::tail(fit$summary.random[[w$idx]], K)   # last K entries = weekly-lag coefs
  keep <- if (is.null(window)) seq_len(K) else which(w$knots >= window[1] & w$knots <= window[2])
  lagspec <- data.table::data.table(
    lag_day = w$knots, OR = exp(at * b$mean),
    low = exp(at * b$`0.025quant`), high = exp(at * b$`0.975quant`)
  )
  smp  <- INLA::inla.posterior.sample(n_sample, fit, verbose = FALSE)
  brow <- utils::tail(grep(paste0("^", w$idx, ":"), rownames(smp[[1]]$latent)), K)[keep]
  cums <- vapply(smp, function(s) exp(at * sum(s$latent[brow, 1])), numeric(1))
  list(lagspec = lagspec,
       cumulative = c(OR = exp(at * sum(b$mean[keep])),
                      low = unname(stats::quantile(cums, 0.025)),
                      high = unname(stats::quantile(cums, 0.975)),
                      p_gt1 = mean(cums > 1)))
}
# back-compat: first weekly var (anseriformes)
bird_rw2_summary <- function(fit, obj, at = 0.5, n_sample = 2000)
  weekly_lag_summary(fit, obj, obj$weekly_vars[1], at = at, n_sample = n_sample)

# Per-week (lag-specific) OR table for the weekly terms: one row per weekly var, wk0..wk4 columns.
weekly_perweek_table <- function(fit, obj, labels = INLA_LABELS) {
  rows <- lapply(obj$weekly_vars, function(v) {
    s <- weekly_lag_summary(fit, obj, v)$lagspec
    data.table::as.data.table(setNames(
      as.list(c(unname(labels[v]), sprintf("%.2f (%.2f, %.2f)", s$OR, s$low, s$high))),
      c("Predictor", paste0("wk", s$lag_day / 7))))
  })
  data.table::rbindlist(rows)
}

# Cumulative-window OR table for the weekly terms: one row per weekly var, one column per window.
weekly_window_table <- function(fit, obj, labels = INLA_LABELS,
                                windows = list("0-28 d" = c(0, 28), "0-14 d" = c(0, 14),
                                               "14-28 d" = c(14, 28))) {
  rows <- lapply(obj$weekly_vars, function(v) {
    cells <- vapply(windows, function(w) {
      cm <- weekly_lag_summary(fit, obj, v, window = w)$cumulative
      sprintf("%.2f (%.2f, %.2f)", cm["OR"], cm["low"], cm["high"])
    }, character(1))
    data.table::as.data.table(setNames(as.list(c(unname(labels[v]), cells)),
                                       c("Predictor", names(windows))))
  })
  data.table::rbindlist(rows)
}

#10. Cumulative OR (per +0.5 SD, or +0.5 SD past the threshold) for every exposure -> one
#    table. at/cen left NULL so each variable gets its own contrast via contrast_for().
cumulative_or_table <- function(fit, obj, at = NULL, cen = NULL, labels = NULL, n_sample = 2000) {
  met_rows <- lapply(obj$met_vars, function(v) {
    cp <- crosspred_cc(fit, obj, v, at = at, cen = cen, cumul = TRUE)
    data.table::data.table(Predictor = v,                    # index by position: one `at` value
                           OR = cp$allRRfit[1], low = cp$allRRlow[1], high = cp$allRRhigh[1])
  })
  # weekly terms are always linear, so they keep the plain +0.5 SD step regardless of
  # whatever threshold contrast the met terms are using
  wk_rows <- lapply(obj$weekly_vars, function(v) {
    cm <- weekly_lag_summary(fit, obj, v, at = 0.5, n_sample = n_sample)$cumulative
    data.table::data.table(Predictor = v, OR = cm["OR"], low = cm["low"], high = cm["high"])
  })
  tbl <- data.table::rbindlist(c(met_rows, wk_rows))
  if (!is.null(labels)) tbl[, Predictor := labels[Predictor]]
  tbl[, `OR (95% CrI)` := sprintf("%.2f (%.2f, %.2f)", OR, low, high)][]
}

#10b. Weekly lag cut-points. Everything gets reported at these four windows, never as one
#     0-28 d cumulative: that lumps a sign-changing lag curve into a single number and
#     inflates the interval for no interpretive gain.
INLA_LAG_WINDOWS <- list("Lag 0-7 days"   = c(0, 7),  "Lag 8-14 days"  = c(8, 14),
                         "Lag 15-21 days" = c(15, 21), "Lag 22-28 days" = c(22, 28))
INLA_LAG_ANCHORS <- c(0, 7, 14, 21, 28)

#10c. OR for one met exposure over a set of lags. Build the contrast vector by hand rather
#     than leaning on crosspred, so a window sum and a single lag go through the same path.
#     crossbasis columns run var-major (v1.l1..v1.lL, v2.l1..), hence the t(outer(v, W)).
met_lag_effect <- function(fit, obj, var, lags, step = NULL) {
  cb <- obj$bases[[var]]
  ct <- contrast_for(obj, var, step)
  ob <- function(x, args) do.call(dlnm::onebasis, c(list(x = x), args))
  v  <- as.numeric(ob(ct$at, attr(cb, "argvar")) - ob(ct$cen, attr(cb, "argvar")))
  W  <- colSums(ob(lags, attr(cb, "arglag")), na.rm = TRUE)
  cvec <- as.vector(t(outer(v, W)))
  cv <- dlnm_coef_vcov(fit, obj$cb_cols[[var]])
  lo <- sum(cvec * cv$coef)
  se <- sqrt(as.numeric(t(cvec) %*% cv$vcov %*% cvec))
  # p_gt1 is the posterior probability of direction, which is what the SI tables report in
  # place of a p-value
  c(OR = exp(lo), low = exp(lo - 1.96 * se), high = exp(lo + 1.96 * se),
    p_gt1 = stats::pnorm(lo / se))
}

#10d. The reporting table: every exposure at the four weekly cut-points, either
#     lag-specific (the OR at that lag) or cumulative (summed across that week's window).
var_max_lag <- function(obj, v) obj$max_lag_by_var[[v]] %||% obj$max_lag

lag_window_table <- function(fit, obj, kind = c("cumulative", "lagspecific"),
                             labels = NULL, n_sample = 2000, step = NULL) {
  kind <- match.arg(kind)
  cols <- if (kind == "cumulative") names(INLA_LAG_WINDOWS) else paste0("Lag ", INLA_LAG_ANCHORS, " days")

  met <- lapply(obj$met_vars, function(v) {
    # a window past this variable's own lag range would be pure extrapolation off the end of
    # its basis, so blank it rather than printing a number nothing supports
    ml <- var_max_lag(obj, v)
    cells <- if (kind == "cumulative")
      vapply(INLA_LAG_WINDOWS, function(w) {
        if (w[2] > ml) return(NA_character_)
        e <- met_lag_effect(fit, obj, v, seq(w[1], w[2]), step)
        sprintf("%.2f (%.2f, %.2f)", e["OR"], e["low"], e["high"]) }, character(1))
    else
      vapply(INLA_LAG_ANCHORS, function(a) {
        if (a > ml) return(NA_character_)
        e <- met_lag_effect(fit, obj, v, a, step)
        sprintf("%.2f (%.2f, %.2f)", e["OR"], e["low"], e["high"]) }, character(1))
    data.table::as.data.table(setNames(as.list(c(v, cells)), c("Predictor", cols)))
  })

  # weekly RW terms carry one coefficient per knot, so a window's cumulative is the sum of
  # the knots falling inside it and the lag-specific is just that knot. these have no stored
  # contrast step of their own, so an unset step means the shared +0.5 SD.
  wstep <- step %||% 0.5
  wk <- lapply(obj$weekly_vars, function(v) {
    s <- weekly_lag_summary(fit, obj, v, at = wstep, n_sample = n_sample)
    cells <- if (kind == "cumulative")
      vapply(INLA_LAG_WINDOWS, function(w) {
        kn <- w[2]                       # one knot per window: 7, 14, 21, 28
        cm <- weekly_lag_summary(fit, obj, v, at = wstep, n_sample = n_sample,
                                 window = c(kn, kn))$cumulative
        sprintf("%.2f (%.2f, %.2f)", cm["OR"], cm["low"], cm["high"]) }, character(1))
    else
      vapply(INLA_LAG_ANCHORS, function(a) {
        r <- s$lagspec[lag_day == a]
        sprintf("%.2f (%.2f, %.2f)", r$OR, r$low, r$high) }, character(1))
    data.table::as.data.table(setNames(as.list(c(v, cells)), c("Predictor", cols)))
  })

  tbl <- data.table::rbindlist(c(met, wk))
  if (!is.null(labels)) tbl[, Predictor := labels[Predictor]]
  tbl[]
}

#10c-bis. Leave-one-out cross-validation. INLA's CPO is the predictive density of each
#     observation with that observation left out, so -sum(log CPO) is the LOO log-score:
#     lower is better, and it is comparable across models fitted to the same data in the same
#     way WAIC is. `failure` flags observations where the LOO approximation is unreliable;
#     a non-zero count means treat the score with care rather than discard it.
loo_score <- function(fit) {
  cpo <- fit$cpo$cpo
  ok  <- is.finite(cpo) & cpo > 0
  list(score = -sum(log(cpo[ok])),
       n_used = sum(ok),
       n_fail = sum(fit$cpo$failure > 0, na.rm = TRUE))
}

#10d-bis. Is a weekly term's lag *shape* actually identified on this panel, or has the RW
#     collapsed to a constant? Measured on the spread of the knot coefficients, which is
#     cleanly bimodal: ~0.004 when flat versus >1.0 when real. A flat one still has a
#     meaningful overall level, so report that once instead of repeating it at every window
#     as though four separate estimates agreed.
weekly_shape_identified <- function(fit, obj, var = obj$weekly_vars[1], tol = 0.05) {
  w <- obj$weekly[[var]]
  if (is.null(w)) return(NA)
  b <- utils::tail(fit$summary.random[[w$idx]], length(w$knots))$mean
  diff(range(b)) > tol
}

#10e. Pr(OR > 1) over the same four weekly windows, same layout as lag_window_table. This is
#     what the SI tables carry instead of a p-value: it's the posterior probability the effect
#     is in the stated direction, so 0.98 and 0.02 are both strong evidence.
prob_window_table <- function(fit, obj, labels = NULL, n_sample = 2000, step = NULL) {
  cols <- names(INLA_LAG_WINDOWS)
  fmt  <- function(p) sprintf("%.2f", p)

  met <- lapply(obj$met_vars, function(v) {
    cells <- vapply(INLA_LAG_WINDOWS, function(w)
      fmt(met_lag_effect(fit, obj, v, seq(w[1], w[2]), step)[["p_gt1"]]), character(1))
    data.table::as.data.table(setNames(as.list(c(v, cells)), c("Predictor", cols)))
  })
  wk <- lapply(obj$weekly_vars, function(v) {
    cells <- vapply(INLA_LAG_WINDOWS, function(w) {
      kn <- w[2]
      fmt(weekly_lag_summary(fit, obj, v, at = step %||% 0.5, n_sample = n_sample,
                             window = c(kn, kn))$cumulative[["p_gt1"]]) }, character(1))
    data.table::as.data.table(setNames(as.list(c(v, cells)), c("Predictor", cols)))
  })

  tbl <- data.table::rbindlist(c(met, wk))
  if (!is.null(labels)) tbl[, Predictor := labels[Predictor]]
  tbl[]
}

#11. Expand config latent-block names (for the dlnm_coef_vcov fallback).
make_latent_names <- function(cfg) {
  unlist(Map(function(tag, l) if (l == 1) tag else paste0(tag, ".", seq_len(l)),
             cfg$contents$tag, cfg$contents$length), use.names = FALSE)
}

# Fixed-effect covariance dug out of the stored config, for when INLA doesn't hand back
# lincomb.derived.covariance.matrix. That happens after inla.core.safe hits a segfault and
# reruns: the fit itself is fine but the correlation matrix is gone, and every cumulative
# contrast needs it. Verified to reproduce the stored matrix to 2e-18 when both exist.
#
# Two traps in here, both silent if you get them wrong:
#   1. contents$start indexes the *full* latent field, which begins with a Predictor block of
#      length nrow(data). Qinv drops that block, so shift everything down by the difference.
#   2. Qinv comes back as a dgCMatrix with only one triangle populated, so half the
#      covariances read as exactly 0. Left unsymmetrised that quietly shrinks contrast SEs by
#      a few percent, which is small enough to look plausible.
config_vcov <- function(fit, cb_cols) {
  cfg <- fit$misc$configs
  if (is.null(cfg)) stop("no covariance matrix and no stored config; refit with config = TRUE")
  Qinv <- cfg$config[[1]]$Qinv
  sel  <- cfg$contents$start[match(cb_cols, cfg$contents$tag)] -
          (sum(cfg$contents$length) - nrow(Qinv))
  if (anyNA(sel) || any(sel < 1 | sel > nrow(Qinv)))
    stop("crossbasis columns don't line up with the stored config")
  V <- as.matrix(Qinv[sel, sel])
  V[V == 0] <- t(V)[V == 0]
  V
}

#12. Exposure labels.
INLA_LABELS <- c(runoff = "Surface Runoff", soil_moisture = "Soil Moisture",
                 temperature = "Temperature", precipitation = "Precipitation",
                 wind_speed = "Wind Speed", snow_cover = "Snow Cover",
                 anseriformes = "Anseriformes", predators = "Predators")

#13. Cumulative-OR forest (one model).
forest_cumulative_ggplot <- function(tbl, title = NULL, labels = INLA_LABELS) {
  d <- data.table::copy(data.table::as.data.table(tbl))
  d[, Predictor := factor(Predictor, levels = rev(unique(Predictor)))]
  ggplot2::ggplot(d, ggplot2::aes(x = OR, y = Predictor)) +
    ggplot2::geom_vline(xintercept = 1, linetype = 2, colour = "gray40") +
    ggplot2::geom_errorbarh(ggplot2::aes(xmin = low, xmax = high), height = 0.2) +
    ggplot2::geom_point(size = 2) +
    ggplot2::scale_x_log10(breaks = c(0.001, 0.01, 0.1, 1, 10, 100),
                           labels = c("0.001", "0.01", "0.1", "1", "10", "100")) +
    ggplot2::labs(x = "Cumulative OR per +0.5 SD (95% CrI)", y = NULL, title = title) +
    theme_spark() + ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 0, hjust = 0.5))
}

#13b. Cumulative-OR forest overlaying a named list of models (e.g. uni vs bi, prior gradient).
forest_group_compare_ggplot <- function(cum_list, title = NULL) {
  grp <- data.table::rbindlist(Map(function(tb, nm)
    data.table::data.table(Group = nm, data.table::as.data.table(tb)[, .(Predictor, OR, low, high)]),
    cum_list, names(cum_list)))
  grp[, Group := factor(Group, levels = names(cum_list))]
  grp[, Predictor := factor(Predictor, levels = rev(unique(cum_list[[1]]$Predictor)))]
  pal <- c("#32135C", "#EAB63E", "#2D6FA3", "#C9405A")[seq_along(cum_list)]
  ggplot2::ggplot(grp, ggplot2::aes(x = OR, y = Predictor, colour = Group)) +
    ggplot2::geom_vline(xintercept = 1, linetype = 2, colour = "gray40") +
    ggplot2::geom_errorbarh(ggplot2::aes(xmin = low, xmax = high), height = 0.25,
                            position = ggplot2::position_dodge(width = 0.6)) +
    ggplot2::geom_point(size = 2, position = ggplot2::position_dodge(width = 0.6)) +
    ggplot2::scale_x_log10(breaks = c(0.001, 0.01, 0.1, 1, 10, 100),
                           labels = c("0.001", "0.01", "0.1", "1", "10", "100")) +
    ggplot2::scale_colour_manual(values = stats::setNames(pal, names(cum_list))) +
    ggplot2::labs(x = "Cumulative OR per +0.5 SD (95% CrI)", y = NULL, title = title) +
    theme_spark() + ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 0, hjust = 0.5))
}

#14. Met lag-response slices (PNG, one panel per exposure).
plot_met_lagresponse_png <- function(fit, obj, file, labels = INLA_LABELS,
                                     width = 11, height = 7) {
  grDevices::png(file, width = width, height = height, units = "in", res = 300)
  on.exit(grDevices::dev.off())
  graphics::par(mfrow = c(ceiling(length(obj$met_vars) / 3), 3), mar = c(4.2, 4.2, 2.5, 1),
                oma = c(3, 4, 0, 1), cex.main = 1.2, cex.lab = 1.2, cex.axis = 1.1,
                tcl = -0.3, mgp = c(2.5, 0.7, 0))
  for (v in obj$met_vars) {
    ct <- contrast_for(obj, v)
    cp <- crosspred_cc(fit, obj, v, cumul = FALSE)
    plot(cp, ptype = "slices", var = ct$at, ci = "area", col = "black",
         ci.arg = list(col = grDevices::rgb(0, 0, 0, 0.1), border = NA),
         lwd = 1.3, main = labels[v], xlab = "", ylab = "", las = 1)
    graphics::abline(h = 1, lty = 2, col = "gray40")
  }
  graphics::mtext("Lag (days)", side = 1, outer = TRUE, line = 0, cex = 1.2)
  graphics::mtext("Odds Ratio (95% CrI)", side = 2, outer = TRUE, line = 0, cex = 1.2)
}

#14b. Overall cumulative exposure-response curve (PNG) for a non-linear crosspred.
plot_overall_png <- function(cp, file, xlab = "Exposure (z-score)", main = "",
                             width = 6, height = 5) {
  grDevices::png(file, width = width, height = height, units = "in", res = 300)
  on.exit(grDevices::dev.off())
  graphics::par(mar = c(4.2, 4.2, 2.5, 1), las = 1, family = "sans")
  plot(cp, "overall", ci = "area", col = "black",
       ci.arg = list(col = grDevices::rgb(0, 0, 0, 0.1), border = NA),
       lwd = 1.6, xlab = xlab, ylab = "Cumulative OR (95% CrI)", main = main)
  graphics::abline(h = 1, lty = 2, col = "gray40")
}

#14c. Exposure x lag OR surface (PNG, PRGn contour) for a non-linear crosspred.
plot_contour_png <- function(cp, file, max_lag, xlab = "Exposure (z-score)",
                             main = "", width = 7, height = 6) {
  y <- cp$predvar
  x <- as.numeric(sub("^lag", "", colnames(cp$matRRfit)))
  z <- t(cp$matRRfit)
  pal    <- rev(RColorBrewer::brewer.pal(11, "PRGn"))
  levels <- base::pretty(z, 20)
  cols   <- c(grDevices::colorRampPalette(pal[1:6])(sum(levels < 1)),
              grDevices::colorRampPalette(pal[6:11])(sum(levels > 1)))
  grDevices::png(file, width = width, height = height, units = "in", res = 300)
  on.exit(grDevices::dev.off())
  graphics::filled.contour(x, y, z, xlab = "Lag (days)", ylab = xlab, main = main,
                           col = cols, levels = levels,
                           plot.axes = { graphics::axis(1); graphics::axis(2) })
}

#15. Weekly bird lag-response (ggplot).
bird_rw2_ggplot <- function(bird_summary, title = NULL) {
  d <- data.table::as.data.table(bird_summary$lagspec)
  ggplot2::ggplot(d, ggplot2::aes(x = lag_day, y = OR)) +
    ggplot2::geom_hline(yintercept = 1, linetype = 2, colour = "gray40") +
    ggplot2::geom_ribbon(ggplot2::aes(ymin = low, ymax = high), alpha = 0.15) +
    ggplot2::geom_line(linewidth = 0.7) + ggplot2::geom_point(size = 1.8) +
    ggplot2::scale_x_continuous(breaks = d$lag_day) +
    ggplot2::labs(x = "Lag (days)", y = "OR per +0.5 SD (95% CrI)",
                  title = title %||% "Anseriformes (weekly lag smooth)") +
    theme_spark() + ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 0, hjust = 0.5))
}

`%||%` <- function(a, b) if (is.null(a)) b else a

