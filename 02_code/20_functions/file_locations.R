#1.Raw data filenames
spillover_data_file = "hpai-outbreaks-clean-placeholder-032124.csv"
fishnet_gpkg        = "fishnet_10km_timeseries.gpkg"
land_cover_gpkg     = "fishnet_landcover.gpkg"

#1b.Multi-state flyway bird-abundance case-crossover panel (eBird joined)
study_abundance_casecrossover_file =
  "study_spillover_events_casecrossover_2022_abundance_060326.csv"

#2.Daily panel and timeseries objects
fishnet_daily_rds             = "mn_fishnet_daily_processed.rds"
fishnet_daily_sf_rds          = "mn_fishnet_daily_processed_sf.rds"
timeseries_daily_rds          = "fishnet_timeseries_daily.RDS"
timeseries_daily_modeling_rds = "fishnet_timeseries_daily_modeling.RDS"
independence_screening_rds    = "independence_screening_daily.RDS"

#3.Climate-anomaly panel (weekly)
climate_anomaly_panel_rds = "climate_anomaly_panel.RDS"

#4.Case-crossover dataframes
casecrossover_rds       = "case_crossover_df.RDS"
casecrossover_uni_rds   = "case_crossover_df_uni.RDS"   # MN unidirectional (primary for INLA DLNM c_11)
casecrossover_1km_rds   = "case_crossover_df_1km.RDS"
casecrossover_200m_rds  = "case_crossover_df_200m.RDS"
casecrossover_14lag_rds = "case_crossover_df_14lag.RDS"
casecrossover_42lag_rds = "case_crossover_df_42lag.RDS"
casecrossover_56lag_rds = "case_crossover_df_56lag.RDS"

#5.Sample-size cascade objects
spillover_summary_rds      = "spillover_summary.RDS"
case_crossover_cascade_rds = "case_crossover_cascade.RDS"

#6.Fitted models
primary_model_rds = "casecrossover_model.RDS"
flex_model_rds    = "casecrossover_flexible_model.RDS"
model_1km_rds     = "casecrossover_model_1km.RDS"
spring_model_rds  = "spring_casecrossover_model.RDS"
fall_model_rds    = "fall_casecrossover_model.RDS"

#7.Multi-state flyway bird case-crossover (c_03) — models and result bundle
flyway_anseri_uni_model_rds = "flyway_anseri_unidirectional.RDS"
flyway_anseri_ts_model_rds  = "flyway_anseri_timestrat.RDS"
flyway_bird_results_rds     = "flyway_bird_casecrossover.RDS"

#8.Multi-state flyway weekly climate+bird case-crossover (c_04)
flyway_weekly_cb_uni_model_rds = "flyway_weekly_climate_bird_unidirectional.RDS"
flyway_weekly_cb_ts_model_rds  = "flyway_weekly_climate_bird_timestrat.RDS"
flyway_weekly_cb_results_rds   = "flyway_weekly_climate_bird.RDS"

#9.Multi-state (7-state) gridded daily+weekly pipeline (a_0*_*_multistate)
#9a.Per-zone weekly bird abundance (a_04 output, consumed by a_01)
multistate_zone_bird_file = "multistate-zone-weekly-bird-abundance-061626.csv"
#9b.Joined panels (a_01)
multistate_fishnet_daily_rds  = "multistate_fishnet_daily_processed.rds"
multistate_fishnet_weekly_rds = "multistate_fishnet_weekly_processed.rds"
multistate_fishnet_sf_rds     = "multistate_fishnet_processed_sf.rds"
#9c.Timeseries panels + independence screening (a_02)
multistate_ts_daily_rds          = "multistate_fishnet_timeseries_daily.RDS"
multistate_ts_daily_modeling_rds = "multistate_fishnet_timeseries_daily_modeling.RDS"
multistate_ts_weekly_rds          = "multistate_fishnet_timeseries_weekly.RDS"
multistate_ts_weekly_modeling_rds = "multistate_fishnet_timeseries_weekly_modeling.RDS"
multistate_independence_screening_rds = "multistate_independence_screening.RDS"
multistate_spillover_summary_rds      = "multistate_spillover_summary.RDS"
#9d.Case-crossover dataframes (a_03). Two case definitions: independent (no suffix,
#    primary = link IND+INR minus spatio-temporal 500 m) and all events (_all).
multistate_casecrossover_daily_rds      = "case_crossover_df_multistate_daily.RDS"
multistate_casecrossover_daily_all_rds  = "case_crossover_df_multistate_daily_all.RDS"
multistate_casecrossover_daily_1km_rds  = "case_crossover_df_multistate_daily_1km.RDS"
multistate_casecrossover_daily_200m_rds = "case_crossover_df_multistate_daily_200m.RDS"
multistate_casecrossover_weekly_rds     = "case_crossover_df_multistate_weekly.RDS"
multistate_casecrossover_weekly_all_rds = "case_crossover_df_multistate_weekly_all.RDS"
multistate_casecrossover_cascade_rds    = "multistate_case_crossover_cascade.RDS"
#9e.Fitted models — independent-only daily (c_05) and weekly (c_06) gridded DLNMs
multistate_model_daily_rds        = "casecrossover_multistate_daily.RDS"
multistate_model_weekly_rds       = "casecrossover_multistate_weekly.RDS"
multistate_model_daily_df2sens_rds = "casecrossover_multistate_daily_df2sens.RDS"  # c_05 §8 sensitivity
