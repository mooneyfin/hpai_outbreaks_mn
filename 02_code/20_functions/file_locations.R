# File locations for HPAI case-crossover analysis (May 2026 daily revision)
# Avoids having to repeat paths every time we load something.

# Confidential spillover data (USDA APHIS)
spillover_data_file = "hpai-outbreaks-clean-placeholder-032124.csv"

# Spatial inputs
fishnet_gpkg = "fishnet_10km_timeseries.gpkg"
land_cover_gpkg = "fishnet_landcover.gpkg"

# Key RDS objects (pre-built for reproducibility)
# Daily panel (primary)
fishnet_daily_rds            = "mn_fishnet_daily_processed.rds"
fishnet_daily_sf_rds         = "mn_fishnet_daily_processed_sf.rds"
timeseries_daily_rds         = "fishnet_timeseries_daily.RDS"
timeseries_daily_modeling_rds= "fishnet_timeseries_daily_modeling.RDS"
independence_screening_rds   = "independence_screening_daily.RDS"

# Climate-anomaly (sensitivity, weekly resolution)
climate_anomaly_panel_rds    = "climate_anomaly_panel.RDS"

# Case-crossover
casecrossover_rds            = "case_crossover_df.RDS"           # primary (500 m × 28 d)
casecrossover_1km_rds        = "case_crossover_df_1km.RDS"       # 1 km sensitivity
casecrossover_200m_rds       = "case_crossover_df_200m.RDS"      # 200 m sensitivity
casecrossover_14lag_rds      = "case_crossover_df_14lag.RDS"     # 0–14 d lag
casecrossover_42lag_rds      = "case_crossover_df_42lag.RDS"     # 0–42 d lag
casecrossover_56lag_rds      = "case_crossover_df_56lag.RDS"     # 0–56 d lag (bird-order)

# Sample-size cascade for CONSORT flowchart
spillover_summary_rds        = "spillover_summary.RDS"
case_crossover_cascade_rds   = "case_crossover_cascade.RDS"

# Model objects
primary_model_rds            = "casecrossover_model.RDS"
flex_model_rds               = "casecrossover_flexible_model.RDS"
model_1km_rds                = "casecrossover_model_1km.RDS"
spring_model_rds             = "spring_casecrossover_model.RDS"
fall_model_rds               = "fall_casecrossover_model.RDS"
