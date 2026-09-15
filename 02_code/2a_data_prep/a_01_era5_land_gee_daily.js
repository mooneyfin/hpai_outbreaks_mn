// ============================================================
// ERA5-Land daily extraction — Mississippi Flyway upper states
// MN + IA + WI + SD + MI + OH + IL + IN, Aug 1 2021 – Dec 31 2022
//
// WARM-UP RE-PULL (Sep 2026): the study period still starts Nov 15 2021, but the
// 90-day rolling anomaly in a_08 needs 90 days of history BEFORE the first study
// day, so it was only complete from Feb 13 2022 and cost 8 cases. The seven
// windows from 2021-08a to 2021-11pre close that gap. Everything from 2021-11a
// onwards already exists in Drive — in the Tasks panel run ONLY the new
// era5_*_2021-08*, -09*, -10* and -11pre exports.
// Output: per-state bi-weekly CSVs + fishnet GeoJSON to Drive/hpai_multistate_era5
//
// Incremental re-pull (Jun 2026): Ohio, Illinois + Indiana added to the original
// 5 states so the climate panel covers the full bird study region. OH/IL already
// have bird abundance; IN still needs a separate eBird extraction before it can
// enter the climate+bird models. The MN/IA/WI/SD/MI CSVs already exist in Drive —
// in the Tasks panel run ONLY the new era5_OH_*, era5_IL_*, era5_IN_* exports,
// plus fishnet_8state_geometry (re-run so a_04/c_04 find the new OH/IL/IN cells).
//
// Memory strategy:
//   * Polygon-mean sampling (reduceRegions)
//   * Bi-weekly export windows (≈15 days each)
//   * Per-state splitting → 8 × 27 = 216 small tasks
//   * 19 bands (added snow_depth, snow_cover, skin_temperature,
//     soil_temperature_level_1, volumetric_soil_water_layer_2, snowfall_hourly;
//     dropped leaf_area_index_high_vegetation, sub_surface_runoff_hourly)
//   * snow_depth is the new primary snow variable; snow_density kept for sensitivity
//   * skin_temperature = land-surface temp (different from 2 m air temp; captures
//     whether surface water bodies are frozen)
//   * soil_temperature_level_1 = 0–7 cm soil temp (frozen vs thawed ground)
//   * volumetric_soil_water_layer_2 = 7–28 cm soil moisture (depth sensitivity)
//   * snowfall_hourly = precip decomposed into snow vs rain (vs total_precipitation)
//   * surface_runoff_hourly retained
// ============================================================

// ---- 0. Configuration ----
var STATES       = ['Minnesota', 'Iowa', 'Wisconsin', 'South Dakota', 'Michigan',
                    'Ohio', 'Illinois', 'Indiana'];
var GRID_M       = 10000;
var DRIVE_FOLDER = 'hpai_multistate_era5';
var EPSG         = 'EPSG:5070';   // NAD83 CONUS Albers (equal-area)

// Bi-weekly export windows. The 2021-08a..2021-11pre rows are warm-up only: they
// never carry a case, they exist so the 90-day baseline is complete on the first
// study day. '2021-11a' below starts Nov 15 because that is the study start; the
// new '2021-11pre' covers the first half of that month.
var windows = [
  ['2021-08-01', '2021-08-15', '2021-08a'],
  ['2021-08-15', '2021-09-01', '2021-08b'],
  ['2021-09-01', '2021-09-15', '2021-09a'],
  ['2021-09-15', '2021-10-01', '2021-09b'],
  ['2021-10-01', '2021-10-15', '2021-10a'],
  ['2021-10-15', '2021-11-01', '2021-10b'],
  ['2021-11-01', '2021-11-15', '2021-11pre'],
  ['2021-11-15', '2021-12-01', '2021-11a'],
  ['2021-12-01', '2021-12-15', '2021-12a'],
  ['2021-12-15', '2022-01-01', '2021-12b'],
  ['2022-01-01', '2022-01-15', '2022-01a'],
  ['2022-01-15', '2022-02-01', '2022-01b'],
  ['2022-02-01', '2022-02-15', '2022-02a'],
  ['2022-02-15', '2022-03-01', '2022-02b'],
  ['2022-03-01', '2022-03-15', '2022-03a'],
  ['2022-03-15', '2022-04-01', '2022-03b'],
  ['2022-04-01', '2022-04-15', '2022-04a'],
  ['2022-04-15', '2022-05-01', '2022-04b'],
  ['2022-05-01', '2022-05-15', '2022-05a'],
  ['2022-05-15', '2022-06-01', '2022-05b'],
  ['2022-06-01', '2022-06-15', '2022-06a'],
  ['2022-06-15', '2022-07-01', '2022-06b'],
  ['2022-07-01', '2022-07-15', '2022-07a'],
  ['2022-07-15', '2022-08-01', '2022-07b'],
  ['2022-08-01', '2022-08-15', '2022-08a'],
  ['2022-08-15', '2022-09-01', '2022-08b'],
  ['2022-09-01', '2022-09-15', '2022-09a'],
  ['2022-09-15', '2022-10-01', '2022-09b'],
  ['2022-10-01', '2022-10-15', '2022-10a'],
  ['2022-10-15', '2022-11-01', '2022-10b'],
  ['2022-11-01', '2022-11-15', '2022-11a'],
  ['2022-11-15', '2022-12-01', '2022-11b'],
  ['2022-12-01', '2022-12-15', '2022-12a'],
  ['2022-12-15', '2023-01-01', '2022-12b']
];

// State name → short tag for filenames
var STATE_TAG = {
  'Minnesota':    'MN',
  'Iowa':         'IA',
  'Wisconsin':    'WI',
  'South Dakota': 'SD',
  'Michigan':     'MI',
  'Ohio':         'OH',
  'Illinois':     'IL',
  'Indiana':      'IN'
};

// ---- 1. State boundaries ----
var statesFc = ee.FeatureCollection('TIGER/2018/States')
  .filter(ee.Filter.inList('NAME', STATES));
var aoi = statesFc.geometry();
Map.addLayer(statesFc, {color: 'red'}, 'AOI states');
Map.centerObject(aoi, 5);

// ---- 2. 10 km fishnet POLYGONS, tagged with state ----
var projection = ee.Projection(EPSG).atScale(GRID_M);
var fishnet_all = aoi.coveringGrid(projection, GRID_M);
fishnet_all = fishnet_all.map(function(f) {
  return f.set('zone_id', f.get('system:index'));
});

// Assign each cell to its state via a spatial join (polygon ↔ polygon intersects).
// Cells whose footprint intersects no state (e.g. over Lake Superior) are
// dropped automatically because saveFirst returns no match.
var spatialFilter = ee.Filter.intersects({
  leftField:  '.geo',
  rightField: '.geo',
  maxError:   100
});
var fishnet_with_state = ee.Join.saveFirst({matchKey: 'state_match'})
  .apply(fishnet_all, statesFc, spatialFilter)
  .map(function(f) {
    var stateFeat = ee.Feature(f.get('state_match'));
    return f.set('state', stateFeat.get('NAME')).set('state_match', null);
  });

print('Total fishnet cells (in 8 states):', fishnet_with_state.size());

// Export the combined fishnet polygons (one GeoJSON for the R pipeline).
// 'fishnet_8state' holds all 8 states' cells, filtered by the `state` column
// downstream. a_01_join_data and c_04 read this file by name — re-run this task
// so the new OH/IL/IN cells are included.
Export.table.toDrive({
  collection:     fishnet_with_state,
  description:    'fishnet_8state_geometry',
  folder:         DRIVE_FOLDER,
  fileNamePrefix: 'fishnet_8state',
  fileFormat:     'GeoJSON'
});

// ---- 3. ERA5-Land variables — 15 bands ----
// Snow representation: snow_depth is the primary snow variable (water-equivalent
// depth in m, ≈ direct proxy for standing snowpack). snow_cover (%) retained as
// a sensitivity. snow_density (kg/m³) kept for backward compatibility with the
// original analysis but is no longer the primary.
var MEAN_VARS = [
  'temperature_2m',
  'dewpoint_temperature_2m',
  'skin_temperature',                  // NEW — land-surface temp (frozen-water proxy)
  'soil_temperature_level_1',          // NEW — 0–7 cm soil temp (frozen vs thawed)
  'u_component_of_wind_10m',
  'v_component_of_wind_10m',
  'surface_pressure',
  'snow_depth',                        // NEW — primary snow variable (m water-eq)
  'snow_cover',                        // NEW — % grid cell with snow
  'snow_density',                      // kept for sensitivity comparison
  'volumetric_soil_water_layer_1',
  'volumetric_soil_water_layer_2',     // NEW — 7–28 cm soil moisture (depth sensitivity)
  'leaf_area_index_low_vegetation'
  // dropped: 'leaf_area_index_high_vegetation'
];
var SUM_VARS = [
  'total_precipitation_hourly',
  'snowfall_hourly',                   // NEW — precip decomposed (snow vs rain)
  'surface_runoff_hourly',
  'total_evaporation_hourly',
  'surface_solar_radiation_downwards_hourly',
  'surface_thermal_radiation_downwards_hourly'
  // dropped: 'sub_surface_runoff_hourly'
];

// ---- 4. Daily aggregation ----
function dailyImage(d) {
  d = ee.Date(d);
  var day = ee.ImageCollection('ECMWF/ERA5_LAND/HOURLY')
    .filterDate(d, d.advance(1, 'day'))
    .filterBounds(aoi);
  var means = day.select(MEAN_VARS).mean();
  var sums  = day.select(SUM_VARS).sum();
  return means.addBands(sums)
    .set('date',              d.format('YYYY-MM-dd'))
    .set('system:time_start', d.millis());
}

// ---- 5. Submit one export per (state, bi-weekly window) ----
STATES.forEach(function(stateName) {
  var stateFishnet = fishnet_with_state.filter(ee.Filter.eq('state', stateName));
  var tag          = STATE_TAG[stateName];

  windows.forEach(function(w) {
    var wStart = ee.Date(w[0]);
    var wEnd   = ee.Date(w[1]);
    var label  = w[2];

    var nDays = wEnd.difference(wStart, 'day');
    var days  = ee.List.sequence(0, nDays.subtract(1))
      .map(function(i) { return wStart.advance(i, 'day'); });
    var dailyImgs = ee.ImageCollection(days.map(dailyImage));

    var sampled = dailyImgs.map(function(img) {
      return img.reduceRegions({
        collection: stateFishnet,
        reducer:    ee.Reducer.mean(),
        scale:      GRID_M,
        crs:        EPSG
      }).map(function(f) {
        return f.set('date', img.get('date')).setGeometry(null);
      });
    }).flatten();

    Export.table.toDrive({
      collection:     sampled,
      description:    'era5_' + tag + '_' + label,
      folder:         DRIVE_FOLDER,
      fileNamePrefix: 'era5_' + tag + '_' + label,
      fileFormat:     'CSV',
      selectors:      ['zone_id', 'state', 'date'].concat(MEAN_VARS).concat(SUM_VARS)
    });
  });
});

print('Submitted', STATES.length * windows.length,
      'per-state bi-weekly exports + 1 fishnet export.');
print('Run all from the Tasks panel.');
