// =====================================================
// ET/PET-BASED VEGETATION WATER-STRESS FRAMEWORK
// WEST AFRICA | 2001–2023
// Author: George Owusu Amoah
// Platform: Google Earth Engine
//
// Primary index:
// ET_PET_WSI = 1 - (ET / PET)
//
// Important interpretation:
// ET_PET_WSI represents evapotranspiration limitation relative
// to atmospheric evaporative demand across mixed vegetated
// landscapes. It is NOT the conventional canopy-temperature CWSI.
// =====================================================


// -----------------------------------------------------
// 1. STUDY AREA AND ANALYSIS SETTINGS
// -----------------------------------------------------
var roi = ee.FeatureCollection(
  'projects/ee-george946/assets/west_Africa'
);

var studyArea = roi.geometry();

Map.centerObject(roi, 4);
Map.addLayer(roi, {}, 'West Africa Boundary', false);

// Full annual series for temporal analysis.
var analysisYears = ee.List.sequence(2001, 2023);

// Benchmark years for mapped outputs and publication figures.
var benchmarkYears = [2001, 2010, 2020, 2023];

// Spatial settings.
var statisticsScale = 5000;
var exportScale = 1000;
var rainfallScale = 6000;

// TVDI edge settings.
var ndviBinWidth = 0.01;
var tvdiEdgeScale = 5000;
var minPixelsPerNDVIBin = 50;


// -----------------------------------------------------
// 2. LANDSAT CLOUD, SHADOW, SNOW, AND SATURATION MASK
// -----------------------------------------------------
function maskLandsatL2(image) {
  var qa = image.select('QA_PIXEL');

  var clearMask = qa.bitwiseAnd(1 << 0).eq(0) // Fill
    .and(qa.bitwiseAnd(1 << 1).eq(0))         // Dilated cloud
    .and(qa.bitwiseAnd(1 << 2).eq(0))         // Cirrus
    .and(qa.bitwiseAnd(1 << 3).eq(0))         // Cloud
    .and(qa.bitwiseAnd(1 << 4).eq(0))         // Cloud shadow
    .and(qa.bitwiseAnd(1 << 5).eq(0));        // Snow

  var saturationMask = image.select('QA_RADSAT').eq(0);

  return image
    .updateMask(clearMask)
    .updateMask(saturationMask);
}


// -----------------------------------------------------
// 3. LANDSAT PREPARATION FUNCTIONS
// -----------------------------------------------------
function prepareLandsat57(image) {
  var red = image.select('SR_B3')
    .multiply(0.0000275)
    .add(-0.2)
    .rename('RED');

  var nir = image.select('SR_B4')
    .multiply(0.0000275)
    .add(-0.2)
    .rename('NIR');

  var lst = image.select('ST_B6')
    .multiply(0.00341802)
    .add(149.0)
    .subtract(273.15)
    .rename('LST');

  return ee.Image.cat([red, nir, lst])
    .copyProperties(image, ['system:time_start']);
}


function prepareLandsat89(image) {
  var red = image.select('SR_B4')
    .multiply(0.0000275)
    .add(-0.2)
    .rename('RED');

  var nir = image.select('SR_B5')
    .multiply(0.0000275)
    .add(-0.2)
    .rename('NIR');

  var lst = image.select('ST_B10')
    .multiply(0.00341802)
    .add(149.0)
    .subtract(273.15)
    .rename('LST');

  return ee.Image.cat([red, nir, lst])
    .copyProperties(image, ['system:time_start']);
}


// -----------------------------------------------------
// 4. LANDSAT COLLECTION FOR EACH YEAR
// -----------------------------------------------------
function getLandsatCollection(year) {
  year = ee.Number(year);

  var start = ee.Date.fromYMD(year, 1, 1);
  var end = start.advance(1, 'year');

  function collectionL57(dataset) {
    return ee.ImageCollection(dataset)
      .filterDate(start, end)
      .filterBounds(studyArea)
      .filter(ee.Filter.eq('PROCESSING_LEVEL', 'L2SP'))
      .map(maskLandsatL2)
      .map(prepareLandsat57);
  }

  function collectionL89(dataset) {
    return ee.ImageCollection(dataset)
      .filterDate(start, end)
      .filterBounds(studyArea)
      .filter(ee.Filter.eq('PROCESSING_LEVEL', 'L2SP'))
      .map(maskLandsatL2)
      .map(prepareLandsat89);
  }

  var l5 = collectionL57('LANDSAT/LT05/C02/T1_L2');
  var l7 = collectionL57('LANDSAT/LE07/C02/T1_L2');
  var l8 = collectionL89('LANDSAT/LC08/C02/T1_L2');
  var l9 = collectionL89('LANDSAT/LC09/C02/T1_L2');

  return l5.merge(l7).merge(l8).merge(l9);
}


// -----------------------------------------------------
// 5. ANNUAL LANDSAT LST AND NDVI
// -----------------------------------------------------
function annualLandsat(year) {
  year = ee.Number(year);

  var composite = getLandsatCollection(year)
    .median()
    .clip(studyArea);

  var rawNDVI = composite.normalizedDifference(['NIR', 'RED']);

  var ndvi = rawNDVI
    .rename('NDVI')
    .updateMask(rawNDVI.gte(0))
    .updateMask(rawNDVI.lte(1));

  var lst = composite.select('LST')
    .rename('LST')
    .updateMask(composite.select('LST').gt(10))
    .updateMask(composite.select('LST').lt(60));

  return ee.Image.cat([lst, ndvi])
    .set({
      'year': year,
      'system:time_start': ee.Date.fromYMD(year, 7, 1).millis()
    })
    .clip(studyArea);
}


// -----------------------------------------------------
// 6. ANNUAL MODIS ET AND PET
// -----------------------------------------------------
function annualETPET(year) {
  year = ee.Number(year);

  var start = ee.Date.fromYMD(year, 1, 1);
  var end = start.advance(1, 'year');

  var modis = ee.ImageCollection('MODIS/061/MOD16A2GF')
    .filterDate(start, end)
    .filterBounds(studyArea)
    .map(function(image) {
      var validMask = image.select('ET').gte(0)
        .and(image.select('PET').gt(0));

      return image.updateMask(validMask);
    });

  var et = modis.select('ET')
    .sum()
    .multiply(0.1)
    .rename('ET')
    .clip(studyArea);

  var pet = modis.select('PET')
    .sum()
    .multiply(0.1)
    .rename('PET')
    .clip(studyArea);

  return ee.Image.cat([et, pet])
    .set({
      'year': year,
      'system:time_start': ee.Date.fromYMD(year, 7, 1).millis()
    })
    .clip(studyArea);
}


// -----------------------------------------------------
// 7. ANNUAL CHIRPS RAINFALL
// -----------------------------------------------------
function annualRainfall(year) {
  year = ee.Number(year);

  var start = ee.Date.fromYMD(year, 1, 1);
  var end = start.advance(1, 'year');

  var rainfall = ee.ImageCollection('UCSB-CHG/CHIRPS/DAILY')
    .filterDate(start, end)
    .filterBounds(studyArea)
    .select('precipitation')
    .sum()
    .rename('RAIN')
    .clip(studyArea);

  return rainfall.set({
    'year': year,
    'system:time_start': ee.Date.fromYMD(year, 7, 1).millis()
  });
}

var annualRainfallCollection = ee.ImageCollection.fromImages(
  analysisYears.map(annualRainfall)
);

var rainfallMean = annualRainfallCollection
  .select('RAIN')
  .mean()
  .rename('RAIN_MEAN');

var rainfallSD = annualRainfallCollection
  .select('RAIN')
  .reduce(ee.Reducer.stdDev())
  .rename('RAIN_SD')
  .max(0.001);

var rainfallZCollection = annualRainfallCollection.map(function(image) {
  var rainfallZ = image.select('RAIN')
    .subtract(rainfallMean)
    .divide(rainfallSD)
    .rename('RAIN_Z');

  return rainfallZ.copyProperties(image, [
    'year',
    'system:time_start'
  ]);
});


// -----------------------------------------------------
// 8. TVDI FROM ANNUAL LST-NDVI SPACE
// Uses 5th and 95th percentile LST values within
// 0.01-wide NDVI bins to estimate wet and dry edges.
// -----------------------------------------------------
function deriveTVDI(lst, ndvi) {
  var ndviBin = ndvi
    .divide(ndviBinWidth)
    .floor()
    .toInt16()
    .rename('NDVI_BIN');

  var reducer = ee.Reducer
    .percentile([5, 95], ['wetEdge', 'dryEdge'])
    .combine(ee.Reducer.count(), '', true)
    .group({
      groupField: 1,
      groupName: 'ndviBin'
    });

  var edgeDictionary = lst.addBands(ndviBin).reduceRegion({
    reducer: reducer,
    geometry: studyArea,
    scale: tvdiEdgeScale,
    maxPixels: 1e13,
    bestEffort: true,
    tileScale: 16
  });

  var edgeList = ee.List(edgeDictionary.get('groups'));

  var edgeFeatures = ee.FeatureCollection(
    edgeList.map(function(item) {
      item = ee.Dictionary(item);

      return ee.Feature(null, {
        ndvi: ee.Number(item.get('ndviBin'))
          .multiply(ndviBinWidth),
        wetEdge: item.get('wetEdge'),
        dryEdge: item.get('dryEdge'),
        count: item.get('count')
      });
    })
  )
  .filter(ee.Filter.notNull(['wetEdge', 'dryEdge']))
  .filter(ee.Filter.gte('count', minPixelsPerNDVIBin))
  .filter(ee.Filter.gte('ndvi', 0))
  .filter(ee.Filter.lte('ndvi', 0.90));

  var wetFit = edgeFeatures.reduceColumns({
    reducer: ee.Reducer.linearFit(),
    selectors: ['ndvi', 'wetEdge']
  });

  var dryFit = edgeFeatures.reduceColumns({
    reducer: ee.Reducer.linearFit(),
    selectors: ['ndvi', 'dryEdge']
  });

  var wetSlope = ee.Number(wetFit.get('scale'));
  var wetOffset = ee.Number(wetFit.get('offset'));

  var drySlope = ee.Number(dryFit.get('scale'));
  var dryOffset = ee.Number(dryFit.get('offset'));

  var lstWet = ndvi.multiply(wetSlope)
    .add(wetOffset)
    .rename('LST_WET');

  var lstDry = ndvi.multiply(drySlope)
    .add(dryOffset)
    .rename('LST_DRY');

  var denominator = lstDry.subtract(lstWet);

  return lst.subtract(lstWet)
    .divide(denominator)
    .updateMask(denominator.abs().gt(0.1))
    .clamp(0, 1)
    .rename('TVDI');
}


// -----------------------------------------------------
// 9. ANNUAL HYDRO-CLIMATIC INDICATORS
// -----------------------------------------------------
function annualAll(year) {
  year = ee.Number(year);

  var landsat = annualLandsat(year);
  var etpet = annualETPET(year);

  var lst = landsat.select('LST');
  var ndvi = landsat.select('NDVI');

  var et = etpet.select('ET');
  var pet = etpet.select('PET');

  var etPetRatio = et.divide(pet)
    .rename('ET_PET')
    .updateMask(pet.gt(0))
    .clamp(0, 1);

  // Primary vegetation water-stress index.
  var etPetWSI = ee.Image(1)
    .subtract(etPetRatio)
    .rename('ET_PET_WSI')
    .clamp(0, 1)
    .updateMask(pet.gt(0));

  // Diagnostic NDVI-only stress expression.
  var ndviStress = ee.Image(1)
    .subtract(ndvi)
    .rename('NDVI_STRESS')
    .clamp(0, 1);

  // Secondary vegetation-weighted test index.
  // This is retained only for the rainfall-performance comparison.
  var vegetationWeightedWSI = etPetWSI
    .multiply(ndviStress)
    .rename('VEG_WEIGHTED_WSI')
    .clamp(0, 1);

  var tvdi = deriveTVDI(lst, ndvi);

  return ee.Image.cat([
    lst,
    ndvi,
    et,
    pet,
    etPetRatio,
    etPetWSI,
    ndviStress,
    vegetationWeightedWSI,
    tvdi
  ])
  .set({
    'year': year,
    'system:time_start': ee.Date.fromYMD(year, 7, 1).millis()
  })
  .clip(studyArea);
}

var annualImages = ee.ImageCollection.fromImages(
  analysisYears.map(annualAll)
).sort('year');

print('Annual Hydro-Climatic Image Collection:', annualImages);


// -----------------------------------------------------
// 10. ADD RAINFALL AND RAINFALL ANOMALIES
// -----------------------------------------------------
var annualWithRain = annualImages.map(function(image) {
  var year = ee.Number(image.get('year'));

  var annualRain = ee.Image(
    annualRainfallCollection
      .filter(ee.Filter.eq('year', year))
      .first()
  );

  var annualRainZ = ee.Image(
    rainfallZCollection
      .filter(ee.Filter.eq('year', year))
      .first()
  );

  return image
    .addBands(annualRain)
    .addBands(annualRainZ)
    .copyProperties(image, [
      'year',
      'system:time_start'
    ]);
});


// -----------------------------------------------------
// 11. MAP VISUALISATION PARAMETERS
// -----------------------------------------------------
var lstVis = {
  min: 25,
  max: 45,
  palette: [
    '#2166ac', '#67a9cf', '#d1e5f0',
    '#fee08b', '#fdae61', '#d73027'
  ]
};

var ndviVis = {
  min: 0,
  max: 0.7,
  palette: [
    '#d73027', '#fdae61', '#fee08b',
    '#a6d96a', '#1a9850'
  ]
};

var etVis = {
  min: 0,
  max: 1200,
  palette: [
    '#ffffcc', '#c2e699', '#78c679',
    '#31a354', '#006837'
  ]
};

var wsiVis = {
  min: 0,
  max: 1,
  palette: [
    '#1a9850', '#91cf60', '#fee08b',
    '#fdae61', '#d73027'
  ]
};

var rainfallZVis = {
  min: -2,
  max: 2,
  palette: [
    '#b2182b', '#ef8a62', '#fddbc7',
    '#f7f7f7', '#d1e5f0', '#67a9cf',
    '#2166ac'
  ]
};


// -----------------------------------------------------
// 12. DISPLAY BENCHMARK-YEAR MAPS
// -----------------------------------------------------
benchmarkYears.forEach(function(year) {
  var image = ee.Image(
    annualWithRain
      .filter(ee.Filter.eq('year', year))
      .first()
  );

  Map.addLayer(
    image.select('LST'),
    lstVis,
    'LST ' + year,
    false
  );

  Map.addLayer(
    image.select('NDVI'),
    ndviVis,
    'NDVI ' + year,
    false
  );

  Map.addLayer(
    image.select('ET'),
    etVis,
    'ET ' + year,
    false
  );

  Map.addLayer(
    image.select('ET_PET_WSI'),
    wsiVis,
    'ET/PET WSI ' + year,
    false
  );

  Map.addLayer(
    image.select('TVDI'),
    wsiVis,
    'TVDI ' + year,
    false
  );

  Map.addLayer(
    image.select('RAIN_Z'),
    rainfallZVis,
    'CHIRPS Rainfall Anomaly ' + year,
    false
  );
});


// -----------------------------------------------------
// 13. ANNUAL REGIONAL SUMMARY STATISTICS
// -----------------------------------------------------
function imageStatistics(year) {
  var image = ee.Image(
    annualWithRain
      .filter(ee.Filter.eq('year', year))
      .first()
  );

  var statistics = image.select([
    'LST',
    'NDVI',
    'ET',
    'ET_PET_WSI',
    'TVDI',
    'RAIN'
  ])
  .reduceRegion({
    reducer: ee.Reducer.mean()
      .combine(ee.Reducer.min(), '', true)
      .combine(ee.Reducer.max(), '', true)
      .combine(ee.Reducer.stdDev(), '', true),
    geometry: studyArea,
    scale: statisticsScale,
    maxPixels: 1e13,
    bestEffort: true,
    tileScale: 16
  });

  return ee.Feature(null, statistics)
    .set('year', year);
}

var statisticsTable = ee.FeatureCollection(
  analysisYears.map(imageStatistics)
).sort('year');

print('Annual Regional Statistics:', statisticsTable);


// -----------------------------------------------------
// 14. SEPARATE TEMPORAL CHARTS
// -----------------------------------------------------
function createTimeSeriesChart(property, title, yAxisTitle) {
  return ui.Chart.feature.byFeature({
    features: statisticsTable,
    xProperty: 'year',
    yProperties: [property]
  })
  .setChartType('LineChart')
  .setOptions({
    title: title,
    hAxis: {title: 'Year'},
    vAxis: {title: yAxisTitle},
    lineWidth: 2,
    pointSize: 4,
    legend: {position: 'none'}
  });
}

print(
  createTimeSeriesChart(
    'LST_mean',
    'Annual Mean Land-Surface Temperature',
    'LST (°C)'
  )
);

print(
  createTimeSeriesChart(
    'NDVI_mean',
    'Annual Mean NDVI',
    'NDVI'
  )
);

print(
  createTimeSeriesChart(
    'ET_mean',
    'Annual Mean Evapotranspiration',
    'ET (mm year⁻¹)'
  )
);

print(
  createTimeSeriesChart(
    'ET_PET_WSI_mean',
    'Annual Mean ET/PET-Based Vegetation Water Stress',
    'ET/PET WSI'
  )
);


// -----------------------------------------------------
// 15. AGGREGATE INDICATORS TO CHIRPS GRID
// Required for pixel-wise rainfall evaluation.
// -----------------------------------------------------
var chirpsProjection = ee.ImageCollection(
  'UCSB-CHG/CHIRPS/DAILY'
)
.first()
.select('precipitation')
.projection();

function aggregateToCHIRPSGrid(image) {
  return image
    .reduceResolution({
      reducer: ee.Reducer.mean(),
      maxPixels: 65535
    })
    .reproject({
      crs: chirpsProjection
    });
}


// -----------------------------------------------------
// 16. PREPARE TIME SERIES FOR RAINFALL EVALUATION
// -----------------------------------------------------
var rainfallEvaluationCollection = annualWithRain.map(function(image) {
  var etPetWSI = aggregateToCHIRPSGrid(
    image.select('ET_PET_WSI')
  ).rename('ET_PET_WSI');

  var tvdi = aggregateToCHIRPSGrid(
    image.select('TVDI')
  ).rename('TVDI');

  var ndviStress = aggregateToCHIRPSGrid(
    image.select('NDVI_STRESS')
  ).rename('NDVI_STRESS');

  var vegetationWeightedWSI = aggregateToCHIRPSGrid(
    image.select('VEG_WEIGHTED_WSI')
  ).rename('VEG_WEIGHTED_WSI');

  return ee.Image.cat([
    etPetWSI,
    tvdi,
    ndviStress,
    vegetationWeightedWSI,
    image.select('RAIN_Z')
  ])
  .copyProperties(image, [
    'year',
    'system:time_start'
  ]);
});


// -----------------------------------------------------
// 17. PIXEL-WISE TEMPORAL CORRELATION WITH RAINFALL
// -----------------------------------------------------
function rainfallCorrelation(indexBand) {
  return rainfallEvaluationCollection
    .select([indexBand, 'RAIN_Z'])
    .reduce(ee.Reducer.pearsonsCorrelation())
    .select('correlation')
    .rename(indexBand + '_RAIN_r');
}

var etPetWSIRainCorrelation = rainfallCorrelation('ET_PET_WSI');
var tvdiRainCorrelation = rainfallCorrelation('TVDI');
var ndviStressRainCorrelation = rainfallCorrelation('NDVI_STRESS');
var vegetationWeightedRainCorrelation = rainfallCorrelation(
  'VEG_WEIGHTED_WSI'
);

var correlationStack = ee.Image.cat([
  etPetWSIRainCorrelation,
  tvdiRainCorrelation,
  ndviStressRainCorrelation,
  vegetationWeightedRainCorrelation
]);

var correlationVis = {
  min: -1,
  max: 1,
  palette: [
    '#b2182b', '#ef8a62', '#fddbc7',
    '#f7f7f7', '#d1e5f0', '#67a9cf',
    '#2166ac'
  ]
};

Map.addLayer(
  etPetWSIRainCorrelation,
  correlationVis,
  'ET/PET WSI–Rainfall Correlation',
  false
);

Map.addLayer(
  tvdiRainCorrelation,
  correlationVis,
  'TVDI–Rainfall Correlation',
  false
);


// -----------------------------------------------------
// 18. MEDIAN CORRELATION SUMMARY
// -----------------------------------------------------
function correlationSummary(name, image) {
  var summary = image.reduceRegion({
    reducer: ee.Reducer.median()
      .combine(ee.Reducer.mean(), '', true)
      .combine(ee.Reducer.stdDev(), '', true),
    geometry: studyArea,
    scale: rainfallScale,
    maxPixels: 1e13,
    bestEffort: true,
    tileScale: 16
  });

  return ee.Feature(null, summary)
    .set('Metric', name);
}

var rainfallPerformanceTable = ee.FeatureCollection([
  correlationSummary(
    'ET/PET-based water-stress index',
    etPetWSIRainCorrelation
  ),
  correlationSummary(
    'Temperature-Vegetation Dryness Index',
    tvdiRainCorrelation
  ),
  correlationSummary(
    'NDVI-only stress expression',
    ndviStressRainCorrelation
  ),
  correlationSummary(
    'Vegetation-weighted ET/PET index',
    vegetationWeightedRainCorrelation
  )
]);

print('Rainfall Evaluation Summary:', rainfallPerformanceTable);


// -----------------------------------------------------
// 19. EXPORT ANNUAL STATISTICS
// -----------------------------------------------------
Export.table.toDrive({
  collection: statisticsTable,
  description: 'WestAfrica_Annual_HydroClimatic_Statistics_2001_2023',
  folder: 'GEE_WestAfrica_WSI',
  fileFormat: 'CSV'
});


// -----------------------------------------------------
// 20. EXPORT RAINFALL PERFORMANCE SUMMARY
// -----------------------------------------------------
Export.table.toDrive({
  collection: rainfallPerformanceTable,
  description: 'WestAfrica_Rainfall_Evaluation_Summary',
  folder: 'GEE_WestAfrica_WSI',
  fileFormat: 'CSV'
});


// -----------------------------------------------------
// 21. EXPORT CORRELATION SAMPLE FOR WILCOXON TEST
// Use the exported CSV in R, Python, or SPSS for the
// paired Wilcoxon signed-rank comparison.
// -----------------------------------------------------
var correlationSamples = correlationStack.sample({
  region: studyArea,
  scale: rainfallScale,
  numPixels: 30000,
  seed: 42,
  geometries: false,
  tileScale: 16
});

Export.table.toDrive({
  collection: correlationSamples,
  description: 'WestAfrica_Rainfall_Correlation_Samples',
  folder: 'GEE_WestAfrica_WSI',
  fileFormat: 'CSV'
});


// -----------------------------------------------------
// 22. EXPORT BENCHMARK-YEAR MULTI-BAND STACKS
// -----------------------------------------------------
benchmarkYears.forEach(function(year) {
  var image = ee.Image(
    annualWithRain
      .filter(ee.Filter.eq('year', year))
      .first()
  );

  Export.image.toDrive({
    image: image.select([
      'LST',
      'NDVI',
      'ET',
      'PET',
      'ET_PET',
      'ET_PET_WSI',
      'TVDI',
      'RAIN',
      'RAIN_Z'
    ]),
    description: 'WestAfrica_HydroClimatic_Stack_' + year,
    folder: 'GEE_WestAfrica_WSI',
    region: studyArea,
    scale: exportScale,
    maxPixels: 1e13
  });
});


// -----------------------------------------------------
// 23. EXPORT RAINFALL-CORRELATION MAPS
// -----------------------------------------------------
Export.image.toDrive({
  image: correlationStack,
  description: 'WestAfrica_Rainfall_Correlation_Maps_2001_2023',
  folder: 'GEE_WestAfrica_WSI',
  region: studyArea,
  scale: rainfallScale,
  maxPixels: 1e13
});


// =====================================================
// END OF SCRIPT
// =====================================================
