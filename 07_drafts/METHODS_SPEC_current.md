# Locked methodology for the GeoHealth revision — as built on 15 September 2026

What the code in `02_code/` fits today and what the exhibits in `04_tables/` and `05_figures/` show.
This is the specification the revised Methods section has to describe. The build-plan narrative in
`REVISION_BUILD_PLAN.md` records how each choice was reached; the annotated manuscript and the response
draft in `07_drafts/revision/` (27 August) predate the final design and are stale.

## 1. Design

| Component | Locked choice |
|---|---|
| Unit of analysis | Cell-day on a 10 km fishnet grid (EPSG:5070), seven northern Mississippi Flyway states |
| Case | The cell on the USDA confirmed-detection date of an independent spillover event |
| Referents | Every other day of the same calendar month in the same cell, excluding the 7 days after detection and any day that was itself a detection. No day-of-week matching. About 23 referents per case |
| Why this design | A fixed calendar partition keeps the conditional likelihood unbiased (Janes, Sheppard & Lumley 2005) and avoids the overlap bias of case-relative referent windows; post-detection days are not valid person-time because the flock is depopulated; conclusions are unchanged for exclusions of 0 to 14 days (`b_03`) |
| Negative control | 25-year climatological normal temperature (1997–2021) as a negative-control exposure (Lipsitch et al. 2010). Minnesota 2.82 under the original unidirectional design, 0.86 under this one; all states 2.99 to 1.13 (Table S4) |
| Independence screen | 218 confirmed detections; 42 removed as not independent introductions (39 farm-to-farm by sequencing and epidemiological investigation, 3 unsequenced); 176 eligible events; 175 case cell-days in 161 cells (Figure S1) |
| Panels | Minnesota 85 cases in 73 cells (primary). Secondary: the other six flyway states (90 cases) and all seven states in one model (175 cases, one set of strata, not pooled estimates) |
| Baseline coverage | ERA5-Land runs from 1 August 2021, so every 2022 case has a full 90-day baseline; no case is lost to the warm-up |
| Known gap | South Dakota is on the grid but has no ERA5-Land export, so it is absent from the flyway comparison despite having 2022 events. State as a limitation |

## 2. Exposures

Six exposures enter every model, defined once in `02_code/20_functions/inla_dlnm_helpers.R` (`INLA_PRIMARY_*`):
temperature, precipitation, soil moisture (0–7 cm), wind speed, surface runoff, and Anseriformes relative abundance.

| Item | Locked choice |
|---|---|
| Primary parameterisation | **90-day shock**: mean of the last 7 days (including the index day) minus the mean of the 90 days ending the day before (day-level version of Zou, Carlton & Grover 2026). Precipitation and runoff are converted to mm and log1p-transformed before differencing |
| Secondary parameterisation | Absolute conditions (daily levels), reported in the supplement (Figure S2, S3, Table S3) |
| Scaling | Each lag column z-scored on the analysis panel; contrasts reported per **+0.5 SD** |
| Birds | eBird Status & Trends weekly relative abundance aggregated to the grid cell, log1p, Anseriformes only. Enters as a level, not a shock |
| Lag window | 0–28 days for every exposure, with a natural-spline lag basis (2 df) |
| Exposure–response | Natural spline (2 df) for each meteorological term; linear for Anseriformes |
| Dropped from the submitted set | Snow density/cover, relative humidity, leaf area index, radiation, total bird abundance, sub-surface runoff, predators (predators kept as a specificity check in Figures 3 and S4) |
| Why precipitation and runoff are both in | The stated hypothesis names them and the interaction model uses them; the exposure set has to be one object across every exhibit even though runoff costs WAIC as a main effect |
| Caveat to state | The temperature shock is 86% explained by day-of-year, so the shock parameterisation must not be described as deseasonalised |

## 3. Model

| Item | Locked choice |
|---|---|
| Likelihood | Conditional Poisson in R-INLA: Poisson outcome with a free intercept per matched set (`iid` term, no sum-to-zero constraint, precision fixed at 1e-6). This reproduces conditional logistic regression; agreement with `clogit` within 0.5% on the ratio scale |
| Crossbases | Built with `dlnm::crossbasis`, entered as fixed effects |
| Prior | N(0, 0.5²) on every crossbasis coefficient, chosen on the implied rate-ratio range for a +0.5 SD contrast (0.66–1.51), not by WAIC (Table S6) |
| Approximation | Simplified Laplace, empirical-Bayes integration |
| Effect measure | **Rate ratio** with 95% **credible interval**, from density sampling of cell-days; plus the posterior probability of direction, reported continuously and never thresholded |
| Reporting | Cumulative rate ratio over four disjoint windows (0–7, 8–14, 15–21, 22–28 d) and the lag-specific contribution at days 0, 7, 14, 21, 28. The full 0–28 d cumulative is not reported |
| Time-invariant covariates | Feedlot density, land cover and every other fixed cell characteristic are conditioned out by design and cannot enter as fixed effects; state the cost (between-cell variation is discarded) |
| Software | R 4.5.3 with R-INLA, dlnm and data.table (no renv.lock is in the repository despite the README; pin versions before resubmission); the case-cell map in Python (geopandas) |

## 4. Interaction analysis (post hoc; `c_07`, Figure 4, Figure S3, Table S9)

Meteorological shock over lag 0–14 d, modified by mean Anseriformes abundance over lag 15–28 d.
One pre-stated specification: strata lag basis with a break at day 8 (independent 0–7 and 8–14 d bins, 0–14 d as their sum from the same posterior);
continuous modifier, centred within stratum and standardised; effects at +1 SD and −1 SD of the modifier; the primary's N(0, 0.5²) prior;
4,000 seeded posterior draws. Fitted for both parameterisations and all three panels. Robust to lag basis, modifier form and bird window (`b_04`, `b_05`; objects in the repository, not tabulated).

## 5. Sensitivity analyses in the supplement

| Exhibit | Content | Script |
|---|---|---|
| Table S1 | eBird species mixture (WQS), R. Kowalski | `b_01` |
| Table S2 | Sample by state | `d_08` |
| Table S3 | Primary results, both parameterisations, three panels, lag-specific and P(dir) | `c_01` |
| Table S4 | Referent design ladder with the negative control: unidirectional (−7 to −28 d), symmetric bidirectional (±14 to 28 d), two-month strata, one-month strata with 7-day post exclusion; Minnesota and all states | `c_02`, `c_03`, `c_04` |
| Table S5 | Lag window and lag basis (Minnesota) | `c_05` |
| Table S6 | Prior on the crossbasis coefficients (Minnesota) | `c_05` |
| Table S7 | Spring versus outside spring (Minnesota) | `c_03` |
| Table S8 | Leave-one-stratum-out cross-validation | `c_06` |
| Table S9 | The interaction, by panel | `c_07` |
| Figures S1–S4 | Selection cascade; forest plot all panels; interaction all panels; waterfowl lag–response all panels and 0–56 d | `d_04`–`d_07` |

Kept runnable but not tabulated: weekly exposure resolution (`b_02`), post-exclusion length (`b_03`), interaction robustness (`b_04`, `b_05`).

## 6. What the submitted Methods text says that no longer holds

- Weekly-aggregated data in 0.09-degree cells; now daily cell-days on the 10 km EPSG:5070 grid, stated one way throughout.
- 101 cases and 110 events with up to four pre-event control weeks; now 85 Minnesota cases in 73 cells (175 across the flyway) with time-stratified referents.
- Conditional logistic regression; now conditional Poisson in INLA with priors. Odds ratios and confidence intervals become rate ratios and credible intervals.
- Predictor set with sub-surface runoff, snow density and total bird abundance; now the six exposures above, with the 90-day shock as the primary parameterisation.
- Lag choice justified by AIC against 6- and 8-week models; now WAIC-selected 0–28 d with a prior-predictive prior (Tables S5, S6).
- Sensitivity list (negative-binomial time series with Duchon spline, season-restricted model, 8-week bird model, 1 km independence screen, LASSO, 2/6/8-week lags, flexible crossbases, mixed-effect anomaly models); replaced by the supplement in section 5.
- Feedlot counts, animal units and NLCD land cover as collected covariates; now conditioned out by design and not in the model.
- The WQS paragraph: still Table S1, but the Methods should say what role it plays now that Anseriformes is the single bird term.

## 7. Decisions still open (confirm before writing)

1. **Contrast scale.** The exhibits use +0.5 SD. The build plan (section 13) argued for per-IQR labelled in natural units; it was never implemented. Keep +0.5 SD or recompute.
2. **One primary or two.** The build plan (section 8) argued both parameterisations are primary; the 15 September commit made Minnesota + 90-day shock the primary with absolute conditions in the supplement. The Methods should say which, and the Discussion has to handle soil moisture changing sign between them.
3. **The log-transform check.** The waterfowl result depended on log1p under the superseded specification; the check on the current specification is not in the supplement list. Decide whether to run and tabulate it.
4. **Power analysis for temperature** (build plan section 20). Not in the supplement list; decide whether it is mentioned.
5. **Timing claims.** The build plan says the waterfowl lag distribution is not identified and bans "peaks at 22–28 days"; Figure 3's per-day curve was built to replace the cumulative curve for that reason. Methods and Results wording must match.
6. **WQS species selection.** Confirm whether the 21-species list still defines the Anseriformes term, or whether the eBird order-level aggregation is used directly.
7. **South Dakota.** Limitation only, or re-export ERA5-Land and rebuild before resubmission.

## 8. Agreed on 14 October 2026

- Decisions 1 to 7 above resolved as follows: keep +0.5 SD and add a natural-units column (state the contrast is evaluated from the stratum mean on the standardised scale); one primary, Minnesota + 90-day shock; rerun the log-transform check on the current spec; include the temperature power analysis; justify Anseriformes a priori from reservoir ecology and cite Table S1 as consistent, non-independent evidence; drop "deseasonalised" (season is controlled by the month strata); run `renv::snapshot()` and commit the lock file.
- **Introduction-to-detection lag (interpretive, no refit).** Index day is the confirmation date. Ssematimba et al. 2024 (PLoS One 19:e0310733) estimate a median 6 d from introduction to first positive sample in US 2022 flocks; Hobbelen et al. 2020 (Sci Rep 10:12388) give 6–15 d (chickens) and 10–19 d (ducks) to clinical detection; our presumptive-to-confirmed lag is median 2 d (95th pct 5, max 10). So introduction sits ~8 d before the index (spread ~5–15 d): lag 0–7 d is post-introduction, 8–14 d is the introduction week, 15–28 d precedes it. Methods state this and the toward-the-null effect of index-date noise on lag-specific estimates; the interaction bins are described as post-introduction week / introduction week / fortnight before; Discussion reads any 0–7 d main effect as detection-side. Do NOT shift the case date or rebuild strata around an assumed introduction day.
- **Direct comparison with Oremush et al. 2026** (Zoonoses Public Health 73:439–450): add a weekly precipitation ≥10 mm indicator sensitivity in the 0–7 and 8–14 d windows; state that relative humidity was tested and null; Discussion contrasts their unidirectional one-control design with the negative-control result (2.82 → 0.86) and notes the opposite sign on precipitation (fall rain in the Fraser Valley vs spring snowmelt in Minnesota).
- Optional: re-anchored lag axis (lags 8–36 d from confirmation) as a supplement table; low priority since Figure S4 already covers lags beyond 28 d.
