# File locations for HPAI case-crossover analysis
# ...avoids having to repeat paths every time we load something!

# Confidential spillover data (USDA APHIS)
spillover_data_file = "hpai-outbreaks-clean-placeholder-032124.csv"

# Spatial inputs
fishnet_gpkg = "fishnet_10km_timeseries.gpkg"
land_cover_gpkg = "fishnet_landcover.gpkg"

# Key RDS objects (pre-built for reproducibility)
timeseries_rds = "fishnet_timeseries.RDS"
timeseries_modeling_rds = "fishnet_timeseries_modeling.RDS"
timeseries_clean_rds = "timeseries_clean.rds"
fishnet_processed_rds = "mn_fishnet_processed.rds"
fishnet_processed_sf_rds = "mn_fishnet_processed_sf.rds"
casecrossover_rds = "case_crossover_df.RDS"
casecrossover_full_rds = "full_case_crossover_df.RDS"
casecrossover_1km_rds = "case_crossover_df_1km.RDS"
casecrossover_200m_rds = "case_crossover_df_200m.RDS"
casecrossover_2lag_rds = "case_crossover_df_2lag.RDS"
casecrossover_6lag_rds = "case_crossover_df_6lag.RDS"
casecrossover_8lag_rds = "case_crossover_df_8lag.RDS"

# Model objects
primary_model_rds = "casecrossover_model.RDS"
flex_model_rds = "casecrossover_flexible_model.RDS"
model_1km_rds = "casecrossover_model_1km.RDS"
spring_model_rds = "spring_casecrossover_model.RDS"
fall_model_rds = "fall_casecrossover_model.RDS"
