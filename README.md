# Associations between Meteorological Conditions, Wild Bird Abundance, and Spillover of Highly Pathogenic Avian Influenza into Poultry Farms: A Case-Crossover Study in Minnesota, 2022

**Authors:** F. M. Mooney, R. Kowalski, C. Gennings, J. L. Warren, K. A. Lehman, N. DeFelice, S. Pei

---

## Project Description

Case-crossover study examining how weekly meteorological conditions and wild bird abundance are associated with HPAI spillover into poultry farms across Minnesota during 2022. We use distributed lag nonlinear models (DLNMs) within a conditional logistic regression framework applied to 10 km grids over 52 weeks.

* Paper: *Submitted to GeoHealth*
* Repository: [https://github.com/sparklabnyc/hpai_outbreaks_mn](https://github.com/sparklabnyc/hpai_outbreaks_mn)

---

## 1. Data

Raw and large processed data files are **not included in this GitHub repository**. They are archived on Zenodo (see DOI below). The key analysis-ready dataframes and fitted model objects are included in this repository as RDS files so that downstream analyses can be reproduced without the full data archive.

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.18983407.svg)](https://doi.org/10.5281/zenodo.18983407)

### 1a. Raw

External input datasets (available on Zenodo):

* **ERA5-Land climate reanalysis** — Weekly meteorological variables (temperature, precipitation, wind, etc.) for Minnesota, 1997–2022. Extracted via Google Earth Engine to the 10 km fishnet grid. Source: [Copernicus Climate Data Store](https://cds.climate.copernicus.eu/)
* **eBird Status & Trends** — Weekly species-level bird abundance estimates. Includes `mn-outbreak-weekly-bird-estimate-2022-012926.csv` (aggregated by outbreak) and `mn-outbreak-weekly-bird-estimate-2022-011626.csv` (earlier version). Source: [Cornell Lab of Ornithology](https://science.ebird.org/en/status-and-trends)
* **Feedlot locations** — Registered feedlot facilities in Minnesota. Source: [Minnesota Pollution Control Agency](https://www.pca.state.mn.us/)
* **Land cover** — National Land Cover Database classifications. Source: [USGS EROS](https://www.usgs.gov/centers/eros/science/national-land-cover-database)
* **HPAI spillover events** — Georeferenced farm locations of confirmed HPAI detections (**confidential — not included in this repository or on Zenodo**). Access must be requested directly from [USDA APHIS](https://www.aphis.usda.gov/livestock-poultry-disease/avian/avian-influenza).

### 1b. Intermediate

Processed spatial and tabular files in `data/processed/`:

* `Case_Zones.gpkg`, `Control_Zones.gpkg`, `Negative_Control_Zones.gpkg` — Case-crossover zone geometries
* `Fishnet_Timeseries_10km.gpkg` — Full 10 km grid with temporal data
* `mn_fishnet_processed.gpkg` / `.csv` — Processed fishnet with joined covariates
* `case_crossover_df_ind*.csv` — Case-crossover datasets at various spatial thresholds
* `mn-hpai-cumulative-lag-case-lag3-control-012926.csv` — Cumulative bird abundance by outbreak, case period lag 3 / control lag +4–+7
* `mn-hpai-cumulative-lag-case-lag4-control-012926.csv` — Cumulative bird abundance by outbreak, case period lag 4 / control lag +4–+8
* `wqs-dataset-guide.csv` — Guide to WQS dataset lag structures

### 1c. Support

* `data/raw/spatial_data/` — 10 km fishnet grid geometry, land cover classification
* `data/file_locations/file_locations.R` — Centralized file path definitions

### 1d. Summary (Analysis-Ready Objects)

Pre-built R objects in `data/objects/` — these are tracked on GitHub so the analysis pipeline can be run immediately:

**Dataframes** (`data/objects/dataframes/`):
* `case_crossover_df.RDS` — Primary case-crossover dataset (4-week lag)
* `case_crossover_df_1km.RDS`, `_200m.RDS` — Spatial restriction variants
* `case_crossover_df_2lag.RDS`, `_6lag.RDS`, `_8lag.RDS` — Lag window variants
* `timeseries_clean.rds` — Clean panel timeseries
* `fishnet_timeseries.RDS` — Full historical timeseries panel (Zenodo only, ~194 MB)
* `fishnet_timeseries_modeling.RDS` — Modeling-ready version (Zenodo only, ~194 MB)

**WQS Results** (`data/objects/`):
* `wqs-model-results-012926.xlsx` — WQS mixture model results (bird species weights and significance)

**Fitted Models** (`data/objects/models/`):
* `casecrossover_model.RDS` — Primary conditional logistic regression model
* `casecrossover_flexible_model.RDS` — Flexible spline variant
* `casecrossover_model_1km.RDS` — 1 km spatial restriction
* `spring_casecrossover_model.RDS`, `fall_casecrossover_model.RDS` — Seasonal models

---

## 2. Code

Each script sources `create_folder_structure.R` and then `code/functions/script_initiate.R`, which loads packages, file locations, objects, and functions.

### Data Preparation (`code/data_prep/`)

These scripts require the confidential spillover data and cannot be run by reviewers. Their outputs are provided as RDS files.

* **a_join_data.Rmd** — Joins spillover events, feedlots, bird abundance, meteorological data, and land cover into the 10 km weekly fishnet panel
* **b_create_timeseries_df.Rmd** — Builds historical climate baselines (10/15/20/25-yr), calculates anomaly z-scores, flags non-independent outbreaks
* **c_create_casscrossover_df.Rmd** — Constructs case-crossover dataset with lagged predictors (0–4 weeks), applies spatial independence exclusions
* **independent_case_confrimation.R** — Validates case independence assumptions
* **d_feedlot_bird_abundance.Rmd** — Extracts eBird species-level abundance at feedlot locations using 3 km rasters (21 species: 15 Anseriformes + 6 non-Anseriformes). Produces `all-mn-feedlot-species-abundance-070125.csv`. Author: Rishi Kowalski

### Models (`code/models/`)

* **a_casecrossover_model.rmd** — Primary DLNM case-crossover model (Figures 3, S3, S5; Tables 3)
* **b_lag_sensitivity.Rmd** — Lag window sensitivity (2, 6, 8 lags)
* **c_season_sensitivity.rmd** — Seasonal stratification (Figures 1, S4; Tables S2, S3)
* **d_bird_order_sensitivity.Rmd** — Bird order stratification: Anseriformes vs. predators (Figure 4; Table S4)
* **e_LASSO_model.Rmd** — Conditional logistic L1 variable selection (Figure S6; Table S5)
* **f_1km_restriction_sensitivity.Rmd** — 1 km spatial independence restriction (Table S6)
* **g_climate_anomaly_models.Rmd** — Climate anomaly models with 10-yr and 25-yr baselines

### Model Processing (`code/model_processing/`)

* **a_table_1.Rmd** — Descriptive Table 1 (Table 1)
* **b_map_casecrossover_df.Rmd** — Case and control zone maps
* **c_case_zone_map.ipynb** — Publication-quality maps in Python/geopandas (Figure 2)
* **d_negative_zone_map.R** — Negative control zone visualization

### Data Exploration (`code/data_exploration/`)

Exploratory analyses not part of the main pipeline:

* **mn-ebird-mixture-analysis.Rmd** — eBird species mixture analysis via WQS (Figure S1; Table S1)
* **casecrossover_validation_models.Rmd** — Leave-one-case-out cross-validation
* **climate_anomaly_model_tests.Rmd** — Climate anomaly model exploration
* **relative_humidity_tests.Rmd** — Relative humidity as a predictor
* **season_model_tests.Rmd** — Seasonal model variants

### Other

* **code/data_acquisition/era5_land_gee.js** — Google Earth Engine script for ERA5-Land extraction
* **code/Archive/** — Earlier script versions (not part of pipeline)

---

## 3. Output

### Manuscript Figure/Table Mapping

| Manuscript Item | Script |
|---|---|
| **Figure 1** | `code/models/c_season_sensitivity.rmd` |
| **Figure 2** | `code/model_processing/c_case_zone_map.ipynb` |
| **Figure 3** | `code/models/a_casecrossover_model.rmd` |
| **Figure 4** | `code/models/d_bird_order_sensitivity.Rmd` |
| **Table 1** | `code/model_processing/a_table_1.Rmd` |
| **Table 3** | `code/models/a_casecrossover_model.rmd` |

### Supplementary Materials

| Item | Script |
|---|---|
| **Figure S1** | `code/data_exploration/mn-ebird-mixture-analysis.Rmd` |
| **Figure S1 (case-ctrl bar/box)** | `code/data_exploration/mn-ebird-mixture-analysis.Rmd` — `mn-case-ctrl-bar-*.png`, `mn-case-ctrl-box-*.png` |
| **Figure S2** | Graphical abstract (not produced by code) |
| **Figure S3** | `code/models/a_casecrossover_model.rmd` |
| **Figure S4** | `code/models/c_season_sensitivity.rmd` |
| **Figure S5** | `code/models/a_casecrossover_model.rmd` |
| **Figure S6** | `code/models/e_LASSO_model.Rmd` |
| **Table S1** | `code/data_exploration/mn-ebird-mixture-analysis.Rmd` |
| **Table S2** | `code/models/c_season_sensitivity.rmd` |
| **Table S3** | `code/models/c_season_sensitivity.rmd` |
| **Table S4** | `code/models/d_bird_order_sensitivity.Rmd` |
| **Table S5** | `code/models/e_LASSO_model.Rmd` |
| **Table S6** | `code/models/f_1km_restriction_sensitivity.Rmd` |

---

## 4. Sensitivity Analyses

* **Lag window** — Models with 2, 6, and 8 lag periods (primary uses 4)
* **Seasonal stratification** — Separate spring and fall models
* **Bird order** — Anseriformes vs. predatory birds
* **Spatial independence** — 1 km exclusion radius around cases
* **Variable selection** — Conditional logistic L1 (LASSO) for reduced model
* **Climate anomaly baselines** — 10-year and 25-year reference periods

---

## Directory Structure

```
hpai_outbreaks_mn/
├── create_folder_structure.R
├── .gitignore
├── README.md
├── code/
│   ├── data_prep/          # Data integration and case-crossover construction
│   ├── models/             # Primary model and sensitivity analyses
│   ├── model_processing/   # Tables, maps, and visualizations
│   ├── data_exploration/   # Exploratory analyses (not in main pipeline)
│   ├── data_acquisition/   # GEE scripts for ERA5-Land extraction
│   ├── functions/          # script_initiate.R, functions.R
│   ├── packages/           # packages_to_load.R
│   └── Archive/            # Earlier script versions
├── data/
│   ├── confidential/       # HPAI spillover locations (restricted)
│   ├── raw/                # External input data (on Zenodo)
│   ├── processed/          # Intermediate spatial and tabular files
│   ├── objects/
│   │   ├── dataframes/     # Analysis-ready RDS objects
│   │   └── models/         # Fitted model RDS objects
│   ├── file_locations/     # file_locations.R
│   └── objects/            # objects.R
├── figures/                # Publication-quality figures
├── tables/                 # Manuscript tables (Word .docx)
├── output/
└── reports/
```

---

## How to Run


1. Download the Zenodo data archive and extract into `data/`, preserving directory structure
2. Open `hpai_outbreaks_mn.Rproj` in RStudio
3. Run any script in `code/models/` or `code/model_processing/` — each script sources `create_folder_structure.R` automatically, which sets up all paths and loads packages

The data preparation scripts (`code/data_prep/`) require the confidential spillover data and cannot be run without it. All downstream analysis scripts work with the pre-built RDS objects provided in the repo.

---

## Dependencies

**R 4.5.2** with key packages:

* Modeling: `survival`, `dlnm`, `clogitL1`, `mgcv`
* Spatial: `sf`, `ggspatial`, `leaflet`
* Data: `dplyr`, `data.table`, `tidyr`, `readr`, `here`
* Visualization: `ggplot2`, `patchwork`, `corrplot`, `cowplot`
* Tables: `gtsummary`, `flextable`, `officer`, `table1`

**Python 3.11** with `geopandas`, `matplotlib`, `contextily` (for `c_case_zone_map.ipynb`).

**Google Earth Engine** (JavaScript API) was used to extract ERA5-Land data. See `code/data_acquisition/era5_land_gee.js`.

---

## Data Availability

| Dataset | Included | Source |
|---|---|---|
| ERA5-Land climate reanalysis | Zenodo | [Copernicus CDS](https://cds.climate.copernicus.eu/) |
| eBird Status & Trends | Zenodo | [Cornell Lab](https://science.ebird.org/en/status-and-trends) |
| MN feedlot locations | Zenodo | [MN Geospatial Commons](https://gisdata.mn.gov/) |
| NLCD land cover | Zenodo | [USGS EROS](https://www.usgs.gov/centers/eros) |
| HPAI spillover events | **Restricted** | Request from [USDA APHIS](https://www.aphis.usda.gov/livestock-poultry-disease/avian/avian-influenza) |
| Analysis-ready RDS objects | GitHub | Produced by `code/data_prep/` |

---

## Citation

Mooney, F. M., Kowalski, R., Gennings, C., Warren, J. L., Lehman, K. A., DeFelice, N., & Pei, S. (2026). Associations between Meteorological Conditions, Wild Bird Abundance, and Spillover of Highly Pathogenic Avian Influenza into Poultry Farms: A Case-Crossover Study in Minnesota, 2022.*Pending Submission*

---

## Contact

* Fintan Mooney
* fm2873@cumc.columbia.edu
* Columbia Mailman School of Public Health
