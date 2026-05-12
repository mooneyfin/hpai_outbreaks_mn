// GEE script to pull ERA5-Land DAILY climate data for the MN HPAI study (revision)
// Daily zone-level extraction Nov 15 2021 – Dec 31 2022 for the time-stratified
// case-crossover analysis at day resolution. Weekly export loop is preserved for
// climatology baselines used in the climate-anomaly sensitivity Rmd only.
//
// Daily outputs → mn_daily_era5land_YYYY.csv (zone_id × date rows)
// Weekly outputs (1997–2021) → mn_weekly_era5land_YYYY.csv (climatology baselines)
//
// Outputs (after Drive download) go to 01_data/1a_exposure_data/meteorological_data/
// See a_00_era5_land_rgee.R for the R-native equivalent.

// ──────────────── Load inputs ────────────────
var fishnet = ee.FeatureCollection('projects/ee-feedlot/assets/mn_fishnet_10km_final');
var era5land = ee.ImageCollection("ECMWF/ERA5_LAND/DAILY_AGGR");

// ──────────────── RH and wind helpers ────────────────
function computeRH(tK, tdK) {
  var tC = tK.subtract(273.15);
  var tdC = tdK.subtract(273.15);
  var es = tC.expression('0.6108 * exp((17.27 * T) / (T + 237.3))', {'T': tC});
  var ea = tdC.expression('0.6108 * exp((17.27 * Td) / (Td + 237.3))', {'Td': tdC});
  return ea.divide(es).multiply(100).rename('relative_humidity');
}
function computeWind(u, v) {
  return u.pow(2).add(v.pow(2)).sqrt().rename('wind_speed');
}

// ──────────────── Fill missing values ────────────────
function fillWithKNN(fishnetCollection, variableList, k) {
  return fishnetCollection.map(function(feature) {
    var neighbors = ee.FeatureCollection(fishnetCollection)
      .filterBounds(feature.geometry().centroid())
      .map(function(f) {
        return f.set('dist', f.geometry().distance(feature.geometry().centroid()));
      })
      .sort('dist');

    var filled = ee.Feature(variableList.iterate(function(v, acc) {
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

    return filled;
  });
}

// ──────────────── Build a single-day zone-level feature collection ────────────────
function dailyZonalForDate(date) {
  var day = ee.Date(date);
  var nextDay = day.advance(1, 'day');
  var dayStr = day.format('YYYYMMdd');

  var eraDay = era5land.filterDate(day, nextDay);

  // DAILY_AGGR already aggregates over the day, so .first() pulls the single image.
  // We still apply .mean() defensively in case of multiple images in window.
  var meanVars = eraDay.mean();
  var sumVars = eraDay.select([
    'total_precipitation_sum',
    'runoff_sum',
    'total_evaporation_sum'
  ]).mean();  // DAILY_AGGR sums are already daily totals; .mean() of one image = same value.

  var t = meanVars.select('temperature_2m');
  var td = meanVars.select('dewpoint_temperature_2m');
  var rh = computeRH(t, td);
  var wind = computeWind(
    meanVars.select('u_component_of_wind_10m'),
    meanVars.select('v_component_of_wind_10m')
  );

  var combined = t.rename('temperature')
    .addBands(td.rename('dewpoint'))
    .addBands(rh)
    .addBands(wind)
    .addBands(meanVars.select('volumetric_soil_water_layer_1').rename('soil_moisture'))
    .addBands(meanVars.select('snow_density'))
    .addBands(sumVars.select('total_precipitation_sum').rename('precipitation'))
    .addBands(sumVars.select('runoff_sum'))
    .addBands(sumVars.select('total_evaporation_sum').rename('evapotranspiration'))
    .addBands(meanVars.select('leaf_area_index_high_vegetation').rename('leaf_area_index'))
    .addBands(meanVars.select('surface_net_solar_radiation_sum').rename('shortwave_radiation'))
    .addBands(meanVars.select('surface_thermal_radiation_downwards_sum').rename('longwave_radiation'))
    .addBands(meanVars.select('surface_pressure'));

  var rawStats = fishnet.map(function(cell) {
    var stats = combined.reduceRegion({
      reducer: ee.Reducer.mean(),
      geometry: cell.geometry(),
      scale: 10000,
      maxPixels: 1e9
    });
    return ee.Feature(cell.geometry())
      .set('zone_id', cell.get('zone_id'))
      .set(stats)
      .set('date', dayStr);
  });

  var filledStats = fillWithKNN(
    rawStats,
    ee.List([
      'temperature', 'dewpoint', 'relative_humidity', 'wind_speed',
      'soil_moisture', 'snow_density', 'precipitation', 'runoff_sum',
      'evapotranspiration', 'leaf_area_index', 'shortwave_radiation',
      'longwave_radiation', 'surface_pressure'
    ]),
    4
  );

  return filledStats.map(function(f) { return f.setGeometry(null); });
}

// ──────────────── Daily exports ────────────────
// Primary analysis window: Nov 15 2021 → Dec 31 2022.
// The Nov–Dec 2021 leading window is required because each case-day's lag-28
// exposure window extends ~28 days before the earliest matching control day in
// its month × DOW stratum, which for early-Feb 2022 cases reaches into Dec 2021.

// 2022 daily
var start2022 = ee.Date('2022-01-01');
var end2022 = ee.Date('2022-12-31');
var nDays2022 = end2022.difference(start2022, 'day').add(1);
var dates2022 = ee.List.sequence(0, nDays2022.subtract(1)).map(function(i) {
  return start2022.advance(i, 'day').format('YYYY-MM-dd');
});

var daily2022 = ee.FeatureCollection(
  dates2022.map(function(d) { return dailyZonalForDate(d); })
).flatten();

Export.table.toDrive({
  collection: daily2022,
  description: 'MN_ERA5Land_Daily_2022',
  folder: 'EarthEngine_MN',
  fileNamePrefix: 'mn_daily_era5land_2022',
  fileFormat: 'CSV'
});

// 2021 trailing window (Nov 15 – Dec 31)
var start2021 = ee.Date('2021-11-15');
var end2021 = ee.Date('2021-12-31');
var nDays2021 = end2021.difference(start2021, 'day').add(1);
var dates2021 = ee.List.sequence(0, nDays2021.subtract(1)).map(function(i) {
  return start2021.advance(i, 'day').format('YYYY-MM-dd');
});

var daily2021 = ee.FeatureCollection(
  dates2021.map(function(d) { return dailyZonalForDate(d); })
).flatten();

Export.table.toDrive({
  collection: daily2021,
  description: 'MN_ERA5Land_Daily_2021_NovDec',
  folder: 'EarthEngine_MN',
  fileNamePrefix: 'mn_daily_era5land_2021_novdec',
  fileFormat: 'CSV'
});


// ──────────────── Weekly climatology baselines (1997–2021) ────────────────
// Retained for the climate-anomaly sensitivity analysis only — NOT used in the
// primary daily case-crossover model. Lives in 02_code/2b_data_exploration/ after revision.
var weeksInYear = ee.List.sequence(0, 51);
var base25yr = ee.List.sequence(1997, 2021);
var base10yr = ee.List.sequence(2012, 2021);

function computeWeeklyClimatology(basePeriod, label) {
  return weeksInYear.map(function(weekIdx) {
    var weekIdxInt = ee.Number(weekIdx);

    var images = basePeriod.map(function(y) {
      var year = ee.Number(y);
      var weekStart = ee.Date.fromYMD(year, 1, 1).advance(weekIdxInt, 'week');
      var weekEnd = weekStart.advance(1, 'week');

      var eraWeek = era5land.filterDate(weekStart, weekEnd);
      var meanVars = eraWeek.mean();
      var sumVars = eraWeek.select([
        'total_precipitation_sum', 'runoff_sum', 'total_evaporation_sum'
      ]).sum();

      var t = meanVars.select('temperature_2m');
      var td = meanVars.select('dewpoint_temperature_2m');
      var rh = computeRH(t, td);
      var wind = computeWind(
        meanVars.select('u_component_of_wind_10m'),
        meanVars.select('v_component_of_wind_10m')
      );

      var combined = t.rename('temperature')
        .addBands(td.rename('dewpoint'))
        .addBands(rh)
        .addBands(wind)
        .addBands(meanVars.select('volumetric_soil_water_layer_1').rename('soil_moisture'))
        .addBands(meanVars.select('snow_density'))
        .addBands(sumVars.select('total_precipitation_sum').rename('precipitation'))
        .addBands(sumVars.select('runoff_sum'))
        .addBands(sumVars.select('total_evaporation_sum').rename('evapotranspiration'))
        .addBands(meanVars.select('leaf_area_index_high_vegetation').rename('leaf_area_index'))
        .addBands(meanVars.select('surface_net_solar_radiation_sum').rename('shortwave_radiation'))
        .addBands(meanVars.select('surface_thermal_radiation_downwards_sum').rename('longwave_radiation'))
        .addBands(meanVars.select('surface_pressure'));

      return combined.set('week_idx', weekIdxInt).set('year', year);
    });

    var imgCol = ee.ImageCollection(images);
    var meanImage = imgCol.mean().set('week_idx', weekIdxInt);
    var sdImage = imgCol.reduce(ee.Reducer.stdDev()).set('week_idx', weekIdxInt);

    var zonal = fishnet.map(function(cell) {
      var meanStats = meanImage.reduceRegion({
        reducer: ee.Reducer.mean(),
        geometry: cell.geometry(),
        scale: 10000,
        maxPixels: 1e9
      });

      var sdStats = sdImage.reduceRegion({
        reducer: ee.Reducer.mean(),
        geometry: cell.geometry(),
        scale: 10000,
        maxPixels: 1e9
      });

      var renamedSdDict = ee.Dictionary(
        sdStats.keys().iterate(function(key, acc) {
          key = ee.String(key);
          var newKey = key.replace('_stdDev', '_sd');
          var updated = ee.Dictionary(acc).set(newKey, sdStats.get(key));
          return updated;
        }, ee.Dictionary({}))
      );

      return cell.set(meanStats)
                 .set(renamedSdDict)
                 .set('week_idx', weekIdxInt)
                 .set('clim_period', label);
    });

    return zonal;
  });
}

var weekly25yr = computeWeeklyClimatology(base25yr, '1997_2021');
var weekly10yr = computeWeeklyClimatology(base10yr, '2012_2021');

var clim25 = ee.FeatureCollection(weekly25yr).flatten();
var clim10 = ee.FeatureCollection(weekly10yr).flatten();

Export.table.toDrive({
  collection: clim25,
  description: 'MN_Climatology_Weekly_25yr',
  folder: 'EarthEngine_MN',
  fileNamePrefix: 'mn_climatology_1997_2021',
  fileFormat: 'CSV'
});

Export.table.toDrive({
  collection: clim10,
  description: 'MN_Climatology_Weekly_10yr',
  folder: 'EarthEngine_MN',
  fileNamePrefix: 'mn_climatology_2012_2021',
  fileFormat: 'CSV'
});


// ──────────────── Create 10km fishnet grid (one-time helper, unchanged) ────────────────
var customBounds = ee.Geometry.Rectangle([-98, 43, -89, 50]);
var cellSize = 0.09;
var lonRange = ee.List.sequence(-98, -89 - cellSize, cellSize);
var latRange = ee.List.sequence(43, 50 - cellSize, cellSize);

var fishnetGrid = ee.FeatureCollection(
  lonRange.map(function(lon) {
    return latRange.map(function(lat) {
      var cell = ee.Geometry.Rectangle([
        lon, lat,
        ee.Number(lon).add(cellSize),
        ee.Number(lat).add(cellSize)
      ]);
      return ee.Feature(cell);
    });
  }).flatten()
);

var fishnetList = fishnetGrid.toList(fishnetGrid.size());
var fishnetWithID = ee.FeatureCollection(
  fishnetList.zip(ee.List.sequence(0, fishnetGrid.size().subtract(1))).map(function(pair) {
    var feature = ee.Feature(ee.List(pair).get(0));
    var index = ee.Number(ee.List(pair).get(1));
    return feature.set('zone_id', index);
  })
);

Export.table.toDrive({
  collection: fishnetWithID,
  description: 'Export_10km_Fishnet_Minnesota',
  fileFormat: 'SHP'
});
