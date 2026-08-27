# Associations between meteorological conditions, wild bird abundance, and spillover of highly pathogenic avian influenza into poultry farms: a Bayesian case-crossover study in Minnesota, 2022
F. M. Mooney, R. Kowalski, C. Gennings, J. L. Warren, K. A. Lehman, N. DeFelice, S. Pei. *Under revision for GeoHealth.*

## Project description

Daily case-crossover analysis of HPAI H5N1 spillover into Minnesota poultry premises in 2022, fitted as a **Bayesian distributed lag model in INLA**. The likelihood is a conditional Poisson with a free intercept per matched set, which reproduces conditional logistic regression while allowing penalised smooths and informative priors.

The unit of analysis is a **cell-day** on a 10 km grid. For each spillover event the case is that cell on the confirmed outbreak date, and referents are the same cell on days −28 to −14 and +14 to +28. Because each cell is its own control, all time-invariant characteristics — poultry density, land cover, terrain, biosecurity — are conditioned out. The **symmetric bidirectional referent window** is a deliberate change from the previous time-stratified design: with referents on both sides of the case a linear seasonal trend cancels exactly, whereas a one-sided design confounds seasonal drift with exposure.

Exposures enter over a 0–28 day lag window on a natural spline lag basis with `N(0, 0.25²)` priors. Meteorological terms are linear on the exposure scale except temperature, modelled as a threshold at 10 °C; weekly waterfowl abundance enters as a random-walk smooth over weekly lag knots and is additionally fitted alone over 0–56 days, since migration operates on a longer timescale than weather.

Analyses run on a common Mississippi Flyway grid covering **seven states with exposure data** — Minnesota, Wisconsin, Iowa, Michigan, Indiana, Illinois and Ohio. Minnesota is the primary population (n = 85 events in 73 cells); the remaining flyway states (n = 90) and the pooled population (n = 175) provide external context. Minnesota and the flyway therefore share one grid, one bird aggregation, and one runoff definition.

> **Known gap.** The grid geometry (`fishnet_8state.geojson`) spans eight states, including South Dakota, because the covering grid was generated over an eight-state area of interest. No ERA5-Land export was ever produced for South Dakota, so it carries no exposure data and never enters the analysis panel. South Dakota recorded HPAI events in 2022, so the flyway comparison omits a state with events for want of exposure data rather than by design. Closing this would mean re-running `a_00_era5_land_gee_multistate.js` for South Dakota and rebuilding the panel; the Minnesota primary analysis is unaffected.

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.18983407.svg)](https://doi.org/10.5281/zenodo.18983407)

## Reproducing the analysis

```
create_folder_structure.R                    # folders
02_code/2a_data_prep/a_0*_*_multistate.*     # build the 8-state daily panel
02_code/2a_data_prep/a_06_build_symmetric_panels.R   # case-crossover panels
02_code/2c_models/c_11_casecrossover_inla_dlnm.Rmd   # primary fit + SI tables S2, S3, S5
02_code/2c_models/c_13_si_extra.Rmd                  # SI tables S4, S7
02_code/2c_models/c_14_crossvalidation.Rmd           # held-out-stratum cross-validation
02_code/2d_model_plotting/d_0[05678]_*               # tables and figures
```

Large intermediates (raw GEE exports, panels, processed rasters) are **not versioned** — see `.gitignore`. Every one is rebuilt by the numbered scripts above from the raw inputs.

## 1. Data

1a_exposure_data: daily ERA5-Land meteorological variables for Minnesota (Nov 15 2021 – Dec 31 2022) extracted via Google Earth Engine; weekly 1997–2021 ERA5-Land for climate-anomaly baselines; weekly eBird Status & Trends abundance estimates at outbreak farms and feedlots.

1b_outcome_data: HPAI spillover events for Minnesota in 2022 (confidential — request from USDA APHIS).

1c_supportive_datasets: 10 km fishnet grid geometry, Minnesota feedlot point shapefile (MPCA), NLCD land-cover classification.

## 2. Code

### 2a. Data prep

a_00_era5_land_gee_multistate.js: Google Earth Engine script to extract daily zone-level ERA5-Land for Nov 15 2021 – Dec 31 2022 across the 8-state Mississippi Flyway grid (10 km fishnet, EPSG:5070). Produces the daily panel used by all current models.

a_00_era5_land_gee_mn8_climatology.js: GEE script for the Minnesota subset of the same grid – weekly 1997–2021 and 2012–2021 climatologies plus the 2022 weekly panel, i.e. the denominator and numerator of the climate-anomaly z-scores.

Superseded (removed, recoverable from git history at commit 1921f09): a_00_era5_land_gee.js and a_00_era5_land_rgee.R, the original 13-band Minnesota-only extraction on the legacy 3,334-cell fishnet. Replaced because they lacked snow_cover and defined runoff as surface + sub-surface, which did not match the flyway panel.

a_01_join_data.Rmd: joins spillover events, daily ERA5-Land, linearly-interpolated daily eBird, feedlots, and NLCD land cover into the daily zone × date panel.

a_02_create_timeseries_df.Rmd: pairwise event distance × day-gap flagging at 200 m / 500 m / 1 km × 28 d; writes daily panel + sample-size cascade.

a_03_create_casecrossover_df.Rmd: builds the time-stratified case-crossover dataframe (zone × year × month × day-of-week), attaches lag 0–28 d exposures, exports spatial and lag-window sensitivity variants.

a_04_feedlot_bird_abundance.Rmd: feedlot-level eBird species abundance extraction (Rishi Kowalski).

a_05_independent_case_confirmation.R: case independence checks against raw farm-level outbreak data.

### 2b. Data exploration

b_00_climate_anomaly_setup.Rmd: 10/15/20/25-year weekly climatology baselines and z-scores for the climate-anomaly sensitivity model.

b_01_ebird_mixture_analysis.Rmd: eBird species mixture analysis via weighted quantile sum (WQS) regression.

b_02_casecrossover_validation.Rmd: leave-one-case-out cross-validation.

b_03_climate_anomaly_tests.Rmd, b_04_relative_humidity_tests.Rmd, b_05_season_model_tests.Rmd: exploratory model variants.

### 2c. Models (main analyses)

c_00_casecrossover_daily.Rmd: daily DLNM case-crossover, lag 0–28 d, conditional logistic regression on zone-day strata.

c_01_casecrossover_weekly.Rmd: weekly DLNM case-crossover (primary analysis), lag 0–4 weeks, conditional logistic regression on 1:4 unidirectional case-control strata.

c_02_survival_cox_ag_weekly.Rmd: Andersen–Gill recurrent-event Cox survival model at weekly resolution.

### 2b. Sensitivity analyses (in `02_code/2b_data_exploration/sensitivity_analyses/`)

Numbered s_NN_*.Rmd scripts covering: lag-window robustness, bird order, LASSO selection, spatial restriction, climate anomalies, crossbasis sweep, daily Cox AG, univariate screen, log-knot lag basis, migration season, Bayesian inference, multi-window timing/power, profile-likelihood bootstrap, control-window variants, daily-vs-weekly comparison, new ERA5-variable swap, poultry-density × runoff interaction (daily and weekly), spring-strata subset.

### 2d. Model plotting

d_00_table_1.Rmd: descriptive Table 1.

d_01_map_casecrossover_df.Rmd: case and control zone maps.

d_02_case_zone_map.ipynb: case zone publication map (Python).

d_03_negative_zone_map.R: negative control zone map.

d_04_consort_flowchart.R: CONSORT-style sample-size cascade.

### 20. Functions

packages_to_load.R: package list with auto-install/load.

file_locations.R: filenames for analysis-ready objects.

script_initiate.R: sources the three files above; called from every analysis script.

functions.R: shared figure palettes and theme helpers (theme_spark, theme_spark_map).

## 3. Output

3a_eda_output: intermediate spatial and tabular files (gpkg, csv).

3b_model_output: analysis-ready dataframes (`dataframes/`) and fitted models (`models/`).

## 4. Tables

Manuscript and supplementary tables (.docx).

## 5. Figures

Manuscript and supplementary figures (.png, .pdf).

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
├── 05_figures
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
| NLCD land cover | Zenodo | [USGS EROS](https://www.usgs.gov/centers/eros) |
| HPAI spillover events | Restricted | Request from [USDA APHIS](https://www.aphis.usda.gov/livestock-poultry-disease/avian/avian-influenza) |
| Analysis-ready RDS objects | GitHub | Produced by `02_code/2a_data_prep/` |

## Contact

Fintan Mooney — fm2873@cumc.columbia.edu — Columbia Mailman School of Public Health
