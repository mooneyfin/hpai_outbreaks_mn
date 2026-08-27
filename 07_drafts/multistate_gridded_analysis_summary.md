# Multi-state (7-state) gridded case-crossover — results summary

*Generated overnight on 2026-06-17. Study region: IA, IL, IN, MI, MN, OH, WI (South Dakota excluded — no eBird extraction). HPAI poultry spillover, 2022.*

## What ran

Full pipeline, end to end, all verified clean (one runner: `02_code/run_multistate_pipeline.R`):

1. **Data prep** — `a_04 → a_01 → a_02 → a_03` (multistate): per-zone weekly bird abundance → 7-state fishnet × day panel (daily ERA5 + weekly→daily smoothed birds) → independence screening → time-stratified case-crossover dataframes.
2. **Models** — `c_05_casecrossover_multistate_gridded.Rmd`: the primary conditional-logistic DLNM (mirroring the MN daily model `c_00`), fit to all four dataframes.

Each step builds **two case definitions**:
- **all events** — every HPAI spillover event (218 outbreaks).
- **independent** — independent introductions only: cross-farm/lateral-transfer (link `CSLT`, 39) and unsequenced (`NOSEQ`, 3) events removed, plus 2 spatio-temporally non-independent second events (within 500 m + 28 d). Leaves 176 link-independent → 173 modeled.

## Sample sizes

| Model | Cases | Controls | Strata | AIC |
|---|---|---|---|---|
| Daily · all events | 213 | 695 | 208 | 602.9 |
| Daily · independent | 173 | 585 | 172 | 498.8 |
| Weekly · all events | 206 | 610 | 191 | 550.5 |
| Weekly · independent | 171 | 550 | 168 | 465.9 |

Daily = lag 0–28 d (zone × year × month × day-of-week strata). Weekly = lag 0–4 wk (zone × year × month strata).

## Cumulative odds ratios (per +0.5 SD, whole lag window)

| Predictor | Daily · all | Daily · independent | Weekly · all | Weekly · independent |
|---|---|---|---|---|
| **Anseriformes** | **2.23 (1.42, 3.50)** | **1.92 (1.23, 3.00)** | **2.96 (1.82, 4.84)** | **2.72 (1.66, 4.46)** |
| Runoff | 1.26 (0.58, 2.73) | 1.68 (0.72, 3.94) | 1.12 (0.72, 1.75) | 1.19 (0.73, 1.94) |
| Soil Moisture | 1.02 (0.44, 2.36) | 1.40 (0.59, 3.34) | 0.88 (0.36, 2.16) | 1.01 (0.40, 2.55) |
| Wind Speed | 1.06 (0.49, 2.31) | 1.11 (0.48, 2.58) | 0.93 (0.60, 1.46) | 0.94 (0.59, 1.49) |
| Temperature | 0.86 (0.48, 1.54) | 0.89 (0.47, 1.68) | 0.81 (0.43, 1.51) | 0.82 (0.42, 1.59) |
| Snow Cover | 0.67 (0.44, 1.02) | **0.58 (0.36, 0.92)** | 0.69 (0.46, 1.01) | 0.69 (0.44, 1.08) |
| Precipitation | 0.61 (0.21, 1.79) | 0.37 (0.11, 1.22) | 0.97 (0.55, 1.72) | 0.87 (0.46, 1.63) |

**Bold** = 95% CI excludes 1.

## Takeaways

1. **Waterfowl (anseriformes) abundance is the dominant, robust risk factor** — significant and positive in *all four* models. Each +0.5 SD raises spillover odds ~2× (daily) to ~3× (weekly). This is the headline result and it does not depend on the case definition.
2. **Snow cover trends protective** (OR ≈ 0.58–0.69), reaching significance in the daily independent model. Plausibly a seasonal/exposure-suppression signal.
3. **The other climate exposures are null here** (precipitation, runoff, soil moisture, temperature, wind speed all cross 1). Pooling 7 states appears to attenuate the weaker climate signals that the MN-only analysis picked up — worth noting as a scale/heterogeneity effect.
4. **All vs. independent barely changes the story.** Adding the 42 cross-farm/unsequenced events shifts point estimates slightly (and tightens CIs via more cases) but flips no conclusion. Waterfowl stays strong either way.
5. **Daily vs. weekly agree in direction**; weekly shows a stronger waterfowl OR (longer aggregation window).

## Outputs to look at

- **Figure**: `05_figures/main/figure_multistate_gridded_cumulative_OR_forest.png` (forest, all vs independent, daily/weekly facets)
- **Table**: `04_tables/main/Table_Multistate_Gridded_Cumulative_OR.docx`
- **Models**: `03_output/3b_model_output/models/casecrossover_multistate_{daily,weekly}_{all,independent}.RDS`
- **Dataframes**: `03_output/3b_model_output/dataframes/case_crossover_df_multistate_{daily,weekly}{,_all}.RDS` (+ daily `_1km`/`_200m` spatial sensitivity)
- **Results bundle**: `…/multistate_gridded_casecrossover_results.RDS`

## Decisions worth a second look

- **"Independent" definition** = link `IND`+`INR` **and** spatio-temporal 500 m removal. If you'd rather it mean link-only, that's a one-line change in `a_03`.
- **Predictor set** matches the MN primary (7 exposures). Could add `predators` or `total_bird_abundance`; both are already in the dataframes as lagged, z-scaled columns.
- **Land cover** is still deferred (NLCD pending). It does **not** affect these case-crossover ORs — it's zone-constant and drops out of the conditional-logit likelihood — so these results are complete without it. Feedlots (MN-only) drop out for the same reason.
- This gridded analysis is **separate** from the existing premises-level weekly flyway models (`c_03`/`c_04`), which were not re-run.
