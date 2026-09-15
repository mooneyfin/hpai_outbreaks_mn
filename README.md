# Associations between meteorological conditions, wild bird abundance, and spillover of highly pathogenic avian influenza into poultry farms: a Bayesian case-crossover study in Minnesota, 2022
F. M. Mooney, R. Kowalski, C. Gennings, J. L. Warren, K. A. Lehman, N. DeFelice, S. Pei. *Under revision for GeoHealth.*

## Project description

Daily case-crossover analysis of HPAI H5N1 spillover into poultry premises across the northern Mississippi Flyway in 2022, fitted as a **Bayesian distributed lag model in INLA**. The likelihood is a conditional Poisson with a free intercept per matched set, which reproduces conditional logistic regression while allowing penalised smooths and informative priors.

The unit of analysis is a **cell-day** on a 10 km grid. For each spillover event the case is that cell on the confirmed detection date; referents are the other days of the same calendar month in the same cell, excluding the seven days after detection (**time-stratified, one-month strata, 7-day post-event exclusion**; about 23 referent cell-days per case). Because each cell is its own control, everything time-invariant about a cell — poultry density, land cover, terrain, biosecurity — is conditioned out. The referent design is checked with a negative-control exposure (the 25-year climatological normal temperature, which cannot cause a 2022 outbreak): it returns 2.99 under the unidirectional design of the original submission and 1.13 under this one (Table S2).

**Six exposures enter every model**: temperature, precipitation, soil moisture, wind speed, runoff and Anseriformes abundance, over a 0–28 day lag. Meteorological terms take a natural-spline exposure–response (2 df) and Anseriformes is linear; the lag basis is a natural spline (2 df); crossbasis coefficients carry `N(0, 0.5²)` priors, chosen on the implied rate-ratio range for a +0.5 SD contrast (Table S6). The **90-day shock** (mean of the last 7 days minus the mean of the preceding 90) is the primary parameterisation; **absolute conditions** are the secondary and sit in the supplement. The exposure set is defined once, in `02_code/20_functions/inla_dlnm_helpers.R` (`INLA_PRIMARY_*`), and read from there by every model.

Panel: **175 spillover cases in 161 cells** — Minnesota 85 in 73 (the primary population), the other six northern Mississippi Flyway states 90, and all seven states together. ERA5-Land runs from 1 August 2021 so every 2022 case has a full 90-day baseline.

The primary analysis is Minnesota under the 90-day shock; the other flyway states, all states together, and absolute conditions are reported in the supplement. The main finding is a rising waterfowl effect through lag 28 days, and — in a post hoc interaction analysis (`c_07`, Figure 4, Table S9) — precipitation and runoff shocks in the fortnight before detection that act only in the wake of high waterfowl abundance two to four weeks earlier.

> **Known gap.** The grid geometry (`fishnet_8state.geojson`) spans eight states, including South Dakota, because the covering grid was generated over an eight-state area of interest. No ERA5-Land export was ever produced for South Dakota, so it carries no exposure data and never enters the analysis panel. South Dakota recorded HPAI events in 2022, so the flyway comparison omits a state with events for want of exposure data rather than by design. Closing this would mean re-running `a_01_era5_land_gee_daily.js` for South Dakota and rebuilding the panel; the Minnesota primary analysis is unaffected.

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.18983407.svg)](https://doi.org/10.5281/zenodo.18983407)

## Reproducing the analysis

R packages are pinned with `renv`; the Python side (the case-cell map) uses the conda environment in `environment.yml` (`conda env create -f environment.yml`). Every script sources `create_folder_structure.R` and `02_code/20_functions/script_initiate.R`.

```
02_code/run_multistate_pipeline.R                    # a_03 -> a_04 -> a_05 -> a_06: the daily panel
02_code/2a_data_prep/a_07_build_referent_panels.R    # referent designs
02_code/2a_data_prep/a_08_build_shock_panel.R        # 90-day shocks; the analysis panel
02_code/2c_models/c_01_primary_models.R              # primary fits, both parameterisations, three panels
02_code/2c_models/c_0{2..7}_*.R                      # sensitivities and the interaction, in order
02_code/2d_model_plotting/d_0{1..7}_*                # figures, in order (d_02 is Python)
02_code/2d_model_plotting/d_0{8,9}_*.Rmd             # tables -> 04_tables/{main,supplement}/*.docx + table_specs/
quarto render 02_code/2d_model_plotting/d_10_manuscript_brief.qmd   # every exhibit, inline, with a read of each
```

Large intermediates (raw GEE exports, panels, fitted models, `03_output/**/*.RDS`) are **not versioned** — see `.gitignore`. Every one is rebuilt by the numbered scripts above from the raw inputs. Scripts are numbered in execution order within each folder.

## 1. Data

1a_exposure_data: daily ERA5-Land meteorological variables for the seven flyway states (1 August 2021 – 31 December 2022) extracted via Google Earth Engine; weekly 1997–2021 ERA5-Land climatology for the same grid (used descriptively in Table 1); weekly eBird Status & Trends abundance estimates aggregated to grid cells.

1b_outcome_data: HPAI spillover events, 2022 (confidential — request from USDA APHIS).

1c_supportive_datasets: 10 km fishnet grid geometry, Minnesota feedlot point shapefile (MPCA), state boundaries.

`01_data/archive/` (not versioned) holds inputs the revised analysis no longer reads: the Minnesota-only ERA5-Land and climatology exports, the Minnesota-only eBird aggregates, the NLCD land-cover rasters, and the legacy Minnesota spatial layers. Land cover was dropped from the analysis because the case-crossover design conditions out everything time-invariant about a cell.

## 2. Code

### 2a. Data prep

a_01_era5_land_gee_daily.js: Google Earth Engine script extracting daily zone-level ERA5-Land across the 8-state flyway grid (10 km fishnet, EPSG:5070).

a_02_era5_land_gee_climatology.js: GEE script for the weekly 1997–2021 climatology on the same grid.

a_03_zone_bird_abundance.Rmd → a_04_join_data.Rmd → a_05_create_timeseries_df.Rmd → a_06_create_casecrossover_df.Rmd: eBird aggregation, the daily zone × date panel, event independence screening, and the case-crossover scaffold. Driven in order by `run_multistate_pipeline.R`.

a_07_build_referent_panels.R: the referent designs compared in Table S2 (symmetric bidirectional, two-month strata, one-month strata with and without post-event exclusion).

a_08_build_shock_panel.R: 90-day baselines and shocks; writes the analysis panel `case_crossover_df_shock_post7.RDS`.

### 2b. Data exploration

b_01_ebird_mixture_analysis.Rmd: eBird species mixture via weighted quantile sum regression (R. Kowalski); the source of Table S1.

Sensitivities kept runnable but not tabulated in the supplement: b_02 (weekly exposure resolution), b_03 (post-event exclusion length), b_04 (the interaction under alternative lag bases and modifier forms) and b_05 (the interaction under alternative bird windows, a weather-first-week-only window, and with each water term dropped). b_05 rebuilds each variant from `c_07` by text substitution, so there is one copy of the interaction model.

### 2c. Models

c_01_primary_models.R: the primary fits — 90-day shock and absolute conditions, Minnesota / other flyway states / all states (Figure 2, Figure S2, Table S3).

c_02_negative_control.R, c_03_unidirectional_and_season.R, c_04_design_ladder.R: referent design ladder with the negative control (Table S4) and the season interaction (Table S7), Minnesota and all states.

c_05_lag_and_prior.R: lag structure and prior sensitivities, all panels (Tables S5, S6 show Minnesota).

c_06_loo_crossvalidation.R: leave-one-stratum-out cross-validation (Table S8).

c_07_interaction.R: meteorological shock over lag 0–14 d modified by Anseriformes abundance over lag 15–28 d — one pre-stated specification (strata lag basis, continuous modifier, the primary's prior, seeded posterior sampling). Figure 4, Figure S3, Table S9.

### 2d. Model plotting

d_01 → d_02 (Python) → d_03: epidemic curve, case-cell map, and their composite (Figure 1). d_04: sample cascade (Figure S1); also writes the case-cell lookup the map draws from. d_05: Figure 2 and Figure S2. d_06: Figure 3 and Figure S4 (waterfowl lag–response). d_07: Figure 4 and Figure S3. d_08: Table 1 (Minnesota sample and conditions by season) and Table S2 (sample by state). d_09: Tables S3–S9. Every table is written to Word and to a display spec that d_10 renders from, so the two cannot drift. d_10_manuscript_brief.qmd: every exhibit in submission order with a short read of each.

`02_code/20_functions/table_helpers.R` is the shared HTML/Word table renderer; `inla_dlnm_helpers.R` holds the model helpers and the canonical exposure set.

### Archive

`02_code/Archive/` holds, by original folder, everything superseded during the revision: the Minnesota-only prep lineage, the unidirectional and symmetric-referent models (c_00–c_16), the climatological-anomaly and time-series threads (c_22, c_38–c_45), and the exploratory notebooks. Kept for provenance; nothing in the manuscript depends on them.

### 20. Functions

packages_to_load.R: package list with auto-install/load.

file_locations.R: filenames for analysis-ready objects.

script_initiate.R: sources the three files above; called from every analysis script.

functions.R: shared figure palettes and theme helpers (theme_spark, theme_spark_map).

## 3. Output

3a_eda_output: intermediate spatial and tabular files (gpkg, csv).

3b_model_output: analysis-ready dataframes (`dataframes/`) and fitted models (`models/`).

## 4. Tables

`main/` Table 1; `supplement/` Tables S2–S9 (.docx, named by exhibit); `table_specs/` the display spec each was exported from, which the manuscript brief renders; `archive/` superseded tables.

## 5. Figures

`main/` Figures 1–4; `supplement/` Figures S1–S4 (.png, named by exhibit); `intermediate/` the epidemic curve and case-cell map that compose Figure 1; `archive/` superseded figures.

## Directory structure

```
.
├── create_folder_structure.R
├── README.md
├── 01_data
│   ├── 1a_exposure_data
│   ├── 1b_outcome_data
│   └── 1c_supportive_datasets
├── 02_code
│   ├── 20_functions
│   ├── 2a_data_prep
│   ├── 2b_data_exploration
│   ├── 2c_models
│   └── 2d_model_plotting
├── 03_output
│   ├── 3a_eda_output
│   └── 3b_model_output
│       ├── dataframes
│       └── models
├── 04_tables
│   ├── main
│   ├── supplement
│   ├── table_specs
│   └── archive
├── 05_figures
│   ├── main
│   ├── supplement
│   ├── intermediate
│   └── archive
├── 06_literature
├── 07_drafts
└── reports
```

Run `create_folder_structure.R` first to create any folders missing on first load.

## Data Availability

| Dataset | Location | Source |
|---|---|---|
| ERA5-Land climate reanalysis | Zenodo | [Copernicus CDS](https://cds.climate.copernicus.eu/) |
| eBird Status & Trends | Zenodo | [Cornell Lab](https://science.ebird.org/en/status-and-trends) |
| Minnesota feedlot locations | Zenodo | [MN Geospatial Commons](https://gisdata.mn.gov/) |
| HPAI spillover events | Restricted | Request from [USDA APHIS](https://www.aphis.usda.gov/livestock-poultry-disease/avian/avian-influenza) |
| Analysis-ready RDS objects | GitHub | Produced by `02_code/2a_data_prep/` |

## Contact

Fintan Mooney — fm2873@cumc.columbia.edu — Columbia Mailman School of Public Health
