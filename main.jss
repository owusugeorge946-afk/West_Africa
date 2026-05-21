// =====================================================
// SATELLITE-BASED CROP WATER STRESS FRAMEWORK
// WEST AFRICA | 2001, 2010, 2020, 2023
// Author: George Owusu Amoah
// Platform: Google Earth Engine
// =====================================================

// -------------------------------
// 1. STUDY AREA
// -------------------------------
var roi = ee.FeatureCollection('projects/ee-george946/assets/west_Africa');
Map.centerObject(roi, 4);
Map.addLayer(roi, {}, 'West Africa Boundary', false);

// -------------------------------
// 2. YEARS
// -------------------------------
var years = [2001, 2010, 2020, 2023];

// -------------------------------
// 3. LANDSAT CLOUD MASK
// -------------------------------
function maskLandsatL2(image) {
  var qa = image.select('QA_PIXEL');

  var cloudShadow = qa.bitwiseAnd(1 << 4).eq(0);
  var clouds = qa.bitwiseAnd(1 << 3).eq(0);
  var cirrus = qa.bitwiseAnd(1 << 2).eq(0);

  return image.updateMask(cloudShadow)
              .updateMask(clouds)
              .updateMask(cirrus);
}

// -------------------------------
// 4. LANDSAT COLLECTION FUNCTION
// -------------------------------
function getLandsatCollection(year) {
  var start = ee.Date.fromYMD(year, 1, 1);
  var end = ee.Date.fromYMD(year, 12, 31);

  var l5 = ee.ImageCollection('LANDSAT/LT05/C02/T1_L2')
    .filterDate(start, end)
    .filterBounds(roi)
    .map(maskLandsatL2)
    .map(function(img) {
      var red = img.select('SR_B3').multiply(0.0000275).add(-0.2).rename('RED');
      var nir = img.select('SR_B4').multiply(0.0000275).add(-0.2).rename('NIR');
      var lst = img.select('ST_B6').multiply(0.00341802).add(149.0).subtract(273.15).rename('LST');
      return img.addBands([red, nir, lst]).select(['RED', 'NIR', 'LST']);
    });

  var l7 = ee.ImageCollection('LANDSAT/LE07/C02/T1_L2')
    .filterDate(start, end)
    .filterBounds(roi)
    .map(maskLandsatL2)
    .map(function(img) {
      var red = img.select('SR_B3').multiply(0.0000275).add(-0.2).rename('RED');
      var nir = img.select('SR_B4').multiply(0.0000275).add(-0.2).rename('NIR');
      var lst = img.select('ST_B6').multiply(0.00341802).add(149.0).subtract(273.15).rename('LST');
      return img.addBands([red, nir, lst]).select(['RED', 'NIR', 'LST']);
    });

  var l8 = ee.ImageCollection('LANDSAT/LC08/C02/T1_L2')
    .filterDate(start, end)
    .filterBounds(roi)
    .map(maskLandsatL2)
    .map(function(img) {
      var red = img.select('SR_B4').multiply(0.0000275).add(-0.2).rename('RED');
      var nir = img.select('SR_B5').multiply(0.0000275).add(-0.2).rename('NIR');
      var lst = img.select('ST_B10').multiply(0.00341802).add(149.0).subtract(273.15).rename('LST');
      return img.addBands([red, nir, lst]).select(['RED', 'NIR', 'LST']);
    });

  var l9 = ee.ImageCollection('LANDSAT/LC09/C02/T1_L2')
    .filterDate(start, end)
    .filterBounds(roi)
    .map(maskLandsatL2)
    .map(function(img) {
      var red = img.select('SR_B4').multiply(0.0000275).add(-0.2).rename('RED');
      var nir = img.select('SR_B5').multiply(0.0000275).add(-0.2).rename('NIR');
      var lst = img.select('ST_B10').multiply(0.00341802).add(149.0).subtract(273.15).rename('LST');
      return img.addBands([red, nir, lst]).select(['RED', 'NIR', 'LST']);
    });

  return l5.merge(l7).merge(l8).merge(l9);
}

// -------------------------------
// 5. ANNUAL LST AND NDVI
// -------------------------------
function annualLandsat(year) {
  var col = getLandsatCollection(year);

  var composite = col.median().clip(roi);

  var ndvi = composite.normalizedDifference(['NIR', 'RED'])
    .rename('NDVI')
    .updateMask(composite.normalizedDifference(['NIR', 'RED']).gte(0));

  var lst = composite.select('LST')
    .rename('LST')
    .updateMask(composite.select('LST').gt(10))
    .updateMask(composite.select('LST').lt(60));

  return ee.Image.cat([lst, ndvi])
    .set('year', year)
    .clip(roi);
}

// -------------------------------
// 6. MODIS ET AND PET
// -------------------------------
function annualETPET(year) {
  var start = ee.Date.fromYMD(year, 1, 1);
  var end = ee.Date.fromYMD(year, 12, 31);

  var modis = ee.ImageCollection('MODIS/061/MOD16A2GF')
    .filterDate(start, end)
    .filterBounds(roi);

  var et = modis.select('ET')
    .sum()
    .multiply(0.1)
    .rename('ET')
    .clip(roi);

  var pet = modis.select('PET')
    .sum()
    .multiply(0.1)
    .rename('PET')
    .clip(roi);

  return ee.Image.cat([et, pet])
    .set('year', year)
    .clip(roi);
}

// -------------------------------
// 7. HYBRID CWSI
// CWSI = (1 - ET/PET) * (1 - NDVI)
// -------------------------------
function annualAll(year) {
  var landsat = annualLandsat(year);
  var etpet = annualETPET(year);

  var lst = landsat.select('LST');
  var ndvi = landsat.select('NDVI');
  var et = etpet.select('ET');
  var pet = etpet.select('PET');

  var et_pet = et.divide(pet)
    .rename('ET_PET')
    .clamp(0, 1);

  var cwsi = ee.Image(1)
    .subtract(et_pet)
    .multiply(ee.Image(1).subtract(ndvi))
    .rename('CWSI')
    .clamp(0, 1)
    .updateMask(ndvi.gte(0))
    .updateMask(pet.gt(0));

  return ee.Image.cat([lst, ndvi, et, pet, et_pet, cwsi])
    .set('year', year)
    .clip(roi);
}

// -------------------------------
// 8. CREATE IMAGE COLLECTION
// -------------------------------
var annualImages = ee.ImageCollection(
  years.map(function(y) {
    return annualAll(y);
  })
);

print('Annual Images:', annualImages);

// -------------------------------
// 9. VISUALISATION PARAMETERS
// -------------------------------
var lstVis = {
  min: 25,
  max: 45,
  palette: ['#2166ac', '#67a9cf', '#d1e5f0', '#fee08b', '#fdae61', '#d73027']
};

var ndviVis = {
  min: 0,
  max: 0.6,
  palette: ['#d73027', '#fdae61', '#fee08b', '#a6d96a', '#1a9850']
};

var etVis = {
  min: 0,
  max: 900,
  palette: ['#ffffcc', '#c2e699', '#78c679', '#31a354', '#006837']
};

var cwsiVis = {
  min: 0,
  max: 1,
  palette: ['#1a9850', '#91cf60', '#fee08b', '#fdae61', '#d73027']
};

// -------------------------------
// 10. DISPLAY MAPS
// -------------------------------
years.forEach(function(year) {
  var img = annualImages.filter(ee.Filter.eq('year', year)).first();

  Map.addLayer(img.select('LST'), lstVis, 'LST ' + year, false);
  Map.addLayer(img.select('NDVI'), ndviVis, 'NDVI ' + year, false);
  Map.addLayer(img.select('ET'), etVis, 'ET ' + year, false);
  Map.addLayer(img.select('CWSI'), cwsiVis, 'CWSI ' + year, false);
});

// -------------------------------
// 11. SUMMARY STATISTICS
// -------------------------------
function imageStats(year) {
  var img = annualImages.filter(ee.Filter.eq('year', year)).first();

  var stats = img.select(['LST', 'NDVI', 'ET', 'CWSI'])
    .reduceRegion({
      reducer: ee.Reducer.mean()
        .combine(ee.Reducer.min(), '', true)
        .combine(ee.Reducer.max(), '', true)
        .combine(ee.Reducer.stdDev(), '', true),
      geometry: roi.geometry(),
      scale: 5000,
      maxPixels: 1e13,
      bestEffort: true,
      tileScale: 16
    });

  return ee.Feature(null, stats).set('year', year);
}

var statsTable = ee.FeatureCollection(
  years.map(function(y) {
    return imageStats(y);
  })
);

print('Summary Statistics:', statsTable);

// -------------------------------
// 12. TEMPORAL CHARTS
// -------------------------------
var chart = ui.Chart.feature.byFeature({
  features: statsTable,
  xProperty: 'year',
  yProperties: ['LST_mean', 'NDVI_mean', 'ET_mean', 'CWSI_mean']
})
.setChartType('LineChart')
.setOptions({
  title: 'Temporal Dynamics of Hydro-Climatic Indicators',
  hAxis: {title: 'Year'},
  vAxis: {title: 'Mean Value'},
  lineWidth: 3,
  pointSize: 6,
  legend: {position: 'bottom'}
});

print(chart);

// -------------------------------
// 13. CORRELATION SAMPLE
// -------------------------------
var sampleImage = annualImages.toBands();

var samples = sampleImage.sample({
  region: roi.geometry(),
  scale: 5000,
  numPixels: 5000,
  seed: 42,
  geometries: false,
  tileScale: 16
});

print('Sample data for correlation/regression:', samples.limit(10));

// -------------------------------
// 14. EXPORT SUMMARY TABLE
// -------------------------------
Export.table.toDrive({
  collection: statsTable,
  description: 'WestAfrica_HydroClimatic_Statistics',
  fileFormat: 'CSV'
});

// -------------------------------
// 15. EXPORT SAMPLE DATA
// -------------------------------
Export.table.toDrive({
  collection: samples,
  description: 'WestAfrica_Correlation_Regression_Samples',
  fileFormat: 'CSV'
});

// -------------------------------
// 16. EXPORT MAPS TO GOOGLE DRIVE
// -------------------------------
years.forEach(function(year) {
  var img = annualImages.filter(ee.Filter.eq('year', year)).first();

  Export.image.toDrive({
    image: img.select('LST'),
    description: 'LST_WestAfrica_' + year,
    folder: 'GEE_WestAfrica_CWSI',
    region: roi.geometry(),
    scale: 1000,
    maxPixels: 1e13
  });

  Export.image.toDrive({
    image: img.select('NDVI'),
    description: 'NDVI_WestAfrica_' + year,
    folder: 'GEE_WestAfrica_CWSI',
    region: roi.geometry(),
    scale: 1000,
    maxPixels: 1e13
  });

  Export.image.toDrive({
    image: img.select('ET'),
    description: 'ET_WestAfrica_' + year,
    folder: 'GEE_WestAfrica_CWSI',
    region: roi.geometry(),
    scale: 1000,
    maxPixels: 1e13
  });

  Export.image.toDrive({
    image: img.select('CWSI'),
    description: 'CWSI_WestAfrica_' + year,
    folder: 'GEE_WestAfrica_CWSI',
    region: roi.geometry(),
    scale: 1000,
    maxPixels: 1e13
  });
});

// -------------------------------
// 17. OPTIONAL: EXPORT STACKED IMAGE
// -------------------------------
years.forEach(function(year) {
  var img = annualImages.filter(ee.Filter.eq('year', year)).first();

  Export.image.toDrive({
    image: img.select(['LST', 'NDVI', 'ET', 'PET', 'ET_PET', 'CWSI']),
    description: 'HydroClimatic_Stack_WestAfrica_' + year,
    folder: 'GEE_WestAfrica_CWSI',
    region: roi.geometry(),
    scale: 1000,
    maxPixels: 1e13
  });
});

// =====================================================
// END OF SCRIPT
// =====================================================
