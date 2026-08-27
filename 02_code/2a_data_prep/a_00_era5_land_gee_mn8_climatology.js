// ERA5-Land weekly climatology + 2022 panel for MINNESOTA, on the 8-state grid.
// only edit section 0, everything below runs off it.
//
//   (A) mn8_climatology_<label>.csv   per zone_id x week_idx, mean + sd of each var
//   (B) mn8_weekly_era5land_2022.csv  the 2022 weekly values, same columns
//
// that's the denominator and the numerator of z = (x - mean) / sd in b_00. the
// daily study-period panel already exists from a_00_era5_land_gee_multistate.js,
// we don't re-export it here.
//
// THE GRID
// zone_id has to match the 8-state panel or nothing joins. section 1 rebuilds it
// exactly the way a_00_era5_land_gee_multistate.js does: cover the FULL 8-state
// aoi, take zone_id from system:index, filter to Minnesota only after that.
// narrowing the aoi first gives you a different grid.
//
// checked offline against fishnet_8state.geojson and against
// multistate_fishnet_timeseries_daily.RDS: this recipe gives 2306 Minnesota cells
// and that zone_id set is identical() to the daily panel's. zone_id looks like
// "-4,227", i.e. x,y in the EPSG:5070 grid, so it's an absolute grid position and
// not a counter. which is also why bolting OH/IL/IN onto the aoi never renumbered
// the original cells.
//
// UNITS
// every column is the weekly MEAN of daily panel values, fluxes included, not a
// weekly sum. b_00 doesn't care either way since both halves are weekly, but
// d_00_table_1.Rmd compares a daily absolute against <var>_mean_25yr, and that
// only means something if the climatology sits on a daily scale. handy side
// effect: every column then matches the 8-state daily panel unit for unit.
//
// derived columns get built per day and then averaged, so they're the mean of the
// daily panel's own values rather than a formula applied to weekly means.
// relative_humidity uses the same Magnus constants as a_01_join_data_multistate.Rmd
// (17.625 / 243.04), not the 17.27 / 237.3 from the old MN-only script. mixing
// the two would bias the RH anomaly.
//
// temperature stays in KELVIN with temperature_C alongside. dewpoint,
// skin_temperature and soil_temperature stay Kelvin too, same as the daily panel.


// ============================================================
// 0. CONFIG, edit here
// ============================================================

// 0a. the aoi. this list, this EPSG and this GRID_M are what define zone_id, so
// they have to stay byte-identical to a_00_era5_land_gee_multistate.js. South
// Dakota included: it's in the fishnet even though the daily panel drops it, and
// pulling it out would relabel the MN/SD border cells.
var STATES = ['Minnesota', 'Iowa', 'Wisconsin', 'South Dakota', 'Michigan',
              'Ohio', 'Illinois', 'Indiana'];
var GRID_M = 10000;
var EPSG   = 'EPSG:5070';   // NAD83 CONUS Albers, equal-area

// 0b. which state's cells to export, and how many there should be. if the count
// below doesn't hit this, stop and fix the grid before anything else.
var TARGET_STATE = 'Minnesota';
var EXPECT_CELLS = 2306;

// 0b2. the fishnet. building the grid in GEE means an 8-way polygon-polygon join
// against TIGER state outlines, which are huge, and it blows the interactive
// memory limit before you can even queue a task. so we don't build it here at
// all: mn8_fishnet_shapefile.zip (next to fishnet_8state.geojson in the repo) is
// the same 2306 cells, already checked identical to the daily panel's zone_ids.
// upload it as a GEE asset once, paste the path here, done.
var FISHNET_ASSET = 'projects/ee-feedlot/assets/mn8_fishnet';

// 0c. climatology baselines. 10yr stays because b_03 and s_06 both want _z_10yr,
// and b_00 reads both. 15yr and 20yr went earlier, nothing consumed them.
var BASELINES = [
  {label: '1997_2021', start: 1997, end: 2021},
  {label: '2012_2021', start: 2012, end: 2021}
];

// 0d. years exported raw, i.e. the years you want anomalies FOR.
var PANEL_YEARS = [2022];

// 0e. 52 weeks of 7 days from Jan 1. week_idx runs 0-51 and it's the join key for
// b_00 and d_00_table_1.Rmd. anchoring on Jan 1 instead of ISO weeks is
// deliberate: isoweek() lands week_idx on a different calendar week every year.
var N_WEEKS = 52;

// 0f. weeks per climatology task. 52 = one file per baseline. drop to 13 if a
// task times out, you'll get 4 files per baseline with a _wNN-NN suffix to rbind.
var WEEKS_PER_TASK = 52;

// 0g. drive folder + filename stems. these land in
// 01_data/1a_exposure_data/meteorological_data/ where b_00 goes looking.
var DRIVE_FOLDER = 'EarthEngine_MN';
var CLIM_PREFIX  = 'mn8_climatology_';
var PANEL_PREFIX = 'mn8_weekly_era5land_';

// 0h. bands. left is the DAILY_AGGR band, right is the column name. names match
// a_01_join_data_multistate.Rmd's raw_to_friendly map so the climatology speaks
// the same vocabulary as the daily panel.
// the multistate script reads ERA5_LAND/HOURLY where flux bands end in _hourly;
// DAILY_AGGR pre-aggregates them to daily totals with a _sum suffix. same
// variables, way cheaper over 25 years.
var MEAN_BANDS = {
  'temperature_2m':                'temperature',        // kelvin
  'dewpoint_temperature_2m':       'dewpoint',           // kelvin
  'skin_temperature':              'skin_temperature',   // kelvin
  'soil_temperature_level_1':      'soil_temperature',   // kelvin
  'surface_pressure':              'surface_pressure',
  'snow_depth':                    'snow_depth',
  'snow_cover':                    'snow_cover',
  'snow_density':                  'snow_density',
  'volumetric_soil_water_layer_1': 'soil_moisture',
  'volumetric_soil_water_layer_2': 'soil_moisture_2',
  'leaf_area_index_low_vegetation':'leaf_area_index'     // LOW, like the 8-state panel
};
var FLUX_BANDS = {
  'total_precipitation_sum':                 'precipitation',
  'snowfall_sum':                            'snowfall',
  'surface_runoff_sum':                      'runoff',   // SURFACE, like the 8-state panel
  'total_evaporation_sum':                   'evapotranspiration',
  'surface_solar_radiation_downwards_sum':   'shortwave_radiation',   // downward
  'surface_thermal_radiation_downwards_sum': 'longwave_radiation'
};

// 0i. gap fill, carried over from the retired a_00_era5_land_gee.js. wired in
// but off, and you
// probably want to leave it off. its filterBounds on the cell centroid only ever
// matches the cell itself on a non-overlapping fishnet, so there's nothing to
// average, and ee.Algorithms.If runs both branches eagerly, so it does a
// reduceColumns over the whole collection for every cell x column x week whether
// anything is null or not. the test block counts nulls for free; flip this on if
// it ever finds one.
var KNN_FILL = false;
var KNN_K    = 4;

// 0j. test first. leave true to queue only the 2022 panel, check its zone_ids
// against the daily panel, then set false to queue the climatologies too.
var RUN_TEST     = true;
var TEST_ONLY    = true;
var TEST_MAP_VAR = 'snow_cover';

// ============================================================
// End of config.
// ============================================================


// ---- 1. the grid ----
// read the uploaded asset. zone_id came from coveringGrid over the full 8-state
// aoi, so these ids line up with the daily panel by construction.
var fishnet  = ee.FeatureCollection(FISHNET_ASSET);
var statesFc = ee.FeatureCollection('TIGER/2018/States')
  .filter(ee.Filter.eq('NAME', TARGET_STATE));

Map.addLayer(fishnet, {color: 'red'}, TARGET_STATE + ' cells');
// centre on the state outline, not on fishnet. centerObject evaluates its
// argument right away and there's no reason to make it chew the whole collection.
Map.centerObject(statesFc.geometry(), 6);

// wrong number here and nothing downstream joins
print('>>> ' + TARGET_STATE + ' cells (must be ' + EXPECT_CELLS + '):', fishnet.size());
print('Sample zone_ids (expect "x,y" like -1,227):',
      fishnet.limit(5).aggregate_array('zone_id'));


// ---- 2. window builders, plain client-side date math ----
function iso(d) { return d.toISOString().slice(0, 10); }

// returns [startISO, stamp, week_idx] for each of the 52 weeks
function weeksOf(year) {
  var out = [];
  for (var i = 0; i < N_WEEKS; i++) {
    var d = new Date(Date.UTC(year, 0, 1));
    d.setUTCDate(d.getUTCDate() + 7 * i);
    out.push([iso(d), iso(d).replace(/-/g, ''), i]);
  }
  return out;
}

function yearsBetween(a, b) {
  var out = [];
  for (var y = a; y <= b; y++) { out.push(y); }
  return out;
}

function weekChunks() {
  var out = [];
  for (var s = 0; s < N_WEEKS; s += WEEKS_PER_TASK) {
    out.push([s, Math.min(s + WEEKS_PER_TASK, N_WEEKS)]);
  }
  return out;
}


// ---- 3. unpack the band config ----
function keysOf(o)      { return Object.keys(o); }
function valuesOf(o)    { return Object.keys(o).map(function(k) { return o[k]; }); }
function suffixed(a, s) { return a.map(function(c) { return c + s; }); }

var MEAN_SRC = keysOf(MEAN_BANDS),  MEAN_OUT = valuesOf(MEAN_BANDS);
var FLUX_SRC = keysOf(FLUX_BANDS),  FLUX_OUT = valuesOf(FLUX_BANDS);
var WIND_SRC = ['u_component_of_wind_10m', 'v_component_of_wind_10m'];
var DERIVED  = ['relative_humidity', 'wind_speed', 'temperature_C'];

var MET_COLS   = MEAN_OUT.concat(FLUX_OUT).concat(DERIVED);
var PANEL_COLS = ['zone_id', 'state', 'week_start', 'week_idx', 'lat', 'lon'].concat(MET_COLS);
var CLIM_COLS  = ['zone_id', 'state', 'week_idx', 'lat', 'lon']
  .concat(suffixed(MET_COLS, '_mean')).concat(suffixed(MET_COLS, '_sd'));

var era5land = ee.ImageCollection('ECMWF/ERA5_LAND/DAILY_AGGR');
var chunks   = weekChunks();

print('Panel columns ('       + PANEL_COLS.length + '):', PANEL_COLS);
print('Climatology columns (' + CLIM_COLS.length  + '):', CLIM_COLS);
print('Rows per file:', N_WEEKS + ' weeks x ' + EXPECT_CELLS + ' cells = ' +
      (N_WEEKS * EXPECT_CELLS));


// ---- 4. one day of the daily panel, as an image ----
// same 19 columns the 8-state daily CSVs carry, plus temperature_C. building the
// daily image first, rather than averaging bands and deriving after, is what
// makes relative_humidity and wind_speed averages of the daily panel's own values.
function dailyPanelImage(d) {
  d = ee.Date(d);
  var img = era5land.filterDate(d, d.advance(1, 'day')).mean();

  var means = img.select(MEAN_SRC).rename(MEAN_OUT);
  var flux  = img.select(FLUX_SRC).rename(FLUX_OUT);
  var uv    = img.select(WIND_SRC);

  // ERA5-Land parks snow_density near 100 kg/m3 when there's no snow, and a_01
  // zeroes that in R. doing it here too is idempotent, and it keeps the baseline
  // on the same scale as the panel without leaning on R to do it.
  var dens = means.select('snow_density');
  means = means.addBands(dens.where(dens.lte(100), 0), ['snow_density'], true);

  var tK  = means.select('temperature');
  var tC  = tK.subtract(273.15);
  var tdC = means.select('dewpoint').subtract(273.15);

  // Magnus, same constants as a_01_join_data_multistate.Rmd
  var rh = tdC.expression(
    '100 * exp((17.625 * Td) / (243.04 + Td)) / exp((17.625 * T) / (243.04 + T))',
    {'Td': tdC, 'T': tC}
  ).rename('relative_humidity');

  var wind = uv.select(WIND_SRC[0]).pow(2)
    .add(uv.select(WIND_SRC[1]).pow(2)).sqrt().rename('wind_speed');

  return means.addBands(flux).addBands(rh).addBands(wind)
    .addBands(tC.rename('temperature_C'));
}

// weekly = mean of the 7 daily panel images
function weeklyImage(startISO) {
  var d0   = ee.Date(startISO);
  var days = ee.List.sequence(0, 6).map(function(i) { return d0.advance(i, 'day'); });
  return ee.ImageCollection(days.map(dailyPanelImage)).mean().rename(MET_COLS);
}

// stack the same week across the baseline years, then mean + sd on the raster.
// b_00 used to do this in R after reducing each year to cells, which isn't
// formally the same sd. a cell is 0.09deg and ERA5-Land is 0.1deg, so it's ~1
// pixel per cell and the gap is rounding, but it isn't identical.
function climatologyImage(weekIdx, years) {
  var stack = ee.ImageCollection(years.map(function(y) {
    return weeklyImage(weeksOf(y)[weekIdx][0]);
  }));
  return stack.mean().rename(suffixed(MET_COLS, '_mean'))
    .addBands(stack.reduce(ee.Reducer.stdDev()).rename(suffixed(MET_COLS, '_sd')));
}


// ---- 5. sampling ----
// crs and scale match a_00_era5_land_gee_multistate.js, so a cell's weekly value
// averages the same pixels its daily values came from.
function sampleCells(img, stamps, keepGeom) {
  return img.reduceRegions({
    collection: fishnet,
    reducer:    ee.Reducer.mean(),
    scale:      GRID_M,
    crs:        EPSG
  }).map(function(f) {
    // cells live in EPSG:5070, so the centroid needs transforming or you get
    // Albers metres instead of lon/lat
    var xy = f.geometry().centroid(1).transform('EPSG:4326', 1).coordinates();
    f = f.set(stamps).set('lon', ee.List(xy).get(0)).set('lat', ee.List(xy).get(1));
    return keepGeom ? f : f.setGeometry(null);
  });
}

function panelCells(entry, keepGeom) {
  return sampleCells(weeklyImage(entry[0]),
                     {week_start: entry[1], week_idx: entry[2]}, keepGeom);
}

function climCells(years) {
  return function(entry, keepGeom) {
    return sampleCells(climatologyImage(entry[2], years), {week_idx: entry[2]}, keepGeom);
  };
}


// ---- 6. gap fill (verbatim from the retired a_00_era5_land_gee.js) ----
function fillWithKNN(fc, variableList, k) {
  return fc.map(function(feature) {
    var neighbors = ee.FeatureCollection(fc)
      .filterBounds(feature.geometry().centroid())
      .map(function(f) {
        return f.set('dist', f.geometry().distance(feature.geometry().centroid()));
      })
      .sort('dist');

    return ee.Feature(variableList.iterate(function(v, acc) {
      v = ee.String(v);
      acc = ee.Feature(acc);
      var val = acc.get(v);
      var fillVal = ee.Algorithms.If(
        ee.Algorithms.IsEqual(val, null),
        neighbors.filter(ee.Filter.neq(v, null)).limit(k)
          .reduceColumns({ reducer: ee.Reducer.mean(), selectors: [v] })
          .get('mean'),
        val
      );
      return acc.set(v, fillVal);
    }, feature));
  });
}

function stack(entries, makeCells, fillCols) {
  var fcs = entries.map(function(e) {
    var fc = makeCells(e, false);
    return KNN_FILL ? fillWithKNN(fc, ee.List(fillCols), KNN_K) : fc;
  });
  return ee.FeatureCollection(fcs).flatten();
}


// ---- 7. test block ----
if (RUN_TEST) {
  var probeWeek = weeksOf(PANEL_YEARS[0])[0];
  var probe     = panelCells(probeWeek, false);

  print('=== TEST: week of', probeWeek[0], '===');
  print('Cells sampled (must be ' + EXPECT_CELLS + '):', probe.size());
  print('Columns present:', probe.first().propertyNames().sort());
  print('First 5 cells:', probe.limit(5));
  print(TEST_MAP_VAR + ' summary:', probe.aggregate_stats(TEST_MAP_VAR));
  print('Cells with any null (0 means KNN_FILL is not needed):',
        probe.size().subtract(probe.filter(ee.Filter.notNull(MET_COLS)).size()));

  Map.addLayer(
    ee.Image().float().paint(panelCells(probeWeek, true), TEST_MAP_VAR),
    {min: 0, max: 100, palette: ['brown', 'tan', 'white', 'lightblue']},
    'TEST ' + TEST_MAP_VAR);
}


// ---- 8. exports ----
function queue(collection, name, selectors) {
  Export.table.toDrive({
    collection:     collection,
    description:    name,
    folder:         DRIVE_FOLDER,
    fileNamePrefix: name,
    fileFormat:     'CSV',
    selectors:      selectors
  });
}

// (B) the 2022 numerator. cheap, and it doubles as the zone_id check.
PANEL_YEARS.forEach(function(y) {
  queue(stack(weeksOf(y), panelCells, MET_COLS), PANEL_PREFIX + y, PANEL_COLS);
});

// (A) the baselines. each one reads its whole baseline, so these are the pricey
// tasks. skipped while TEST_ONLY is on.
var nClim = 0;
if (!TEST_ONLY) {
  BASELINES.forEach(function(base) {
    var years = yearsBetween(base.start, base.end);
    chunks.forEach(function(ch) {
      var suffix = chunks.length > 1
        ? '_w' + ('0' + ch[0]).slice(-2) + '-' + ('0' + (ch[1] - 1)).slice(-2)
        : '';
      // weeksOf() is just an index carrier here: climCells reads entry[2] and
      // climatologyImage works out the per-year start dates itself
      queue(stack(weeksOf(base.end).slice(ch[0], ch[1]), climCells(years),
                  suffixed(MET_COLS, '_mean').concat(suffixed(MET_COLS, '_sd'))),
            CLIM_PREFIX + base.label + suffix, CLIM_COLS);
      nClim++;
    });
  });
}

print('Queued', PANEL_YEARS.length + nClim, 'exports.');
print(TEST_ONLY
  ? '>>> TEST_ONLY is on: only ' + PANEL_PREFIX + PANEL_YEARS[0] + '. Check its ' +
    'zone_ids against the daily panel, then set TEST_ONLY = false.'
  : 'Run them from the Tasks panel.');
