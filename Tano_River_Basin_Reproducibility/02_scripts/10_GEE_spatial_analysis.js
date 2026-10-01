// Tano River Basin spatial regeneration workflow.
// CDD = maximum annual consecutive run of daily precipitation < 1 mm.
var basin=ee.FeatureCollection('projects/ee-george946/assets/Tano_River');
var geom=basin.geometry(); var years=[2001,2010,2020,2023]; Map.centerObject(basin,8);
var chirps=ee.ImageCollection('UCSB-CHG/CHIRPS/DAILY').select('precipitation');
function annualCHIRPS(y){
 y=ee.Number(y); var s=ee.Date.fromYMD(y,1,1),e=s.advance(1,'year'); var c=chirps.filterDate(s,e);
 var rain=c.sum().rename('rainfall'), rx1=c.max().rename('rx1day'), list=c.sort('system:time_start').toList(c.size());
 var init=ee.Dictionary({run:ee.Image(0).toInt16(),maxrun:ee.Image(0).toInt16(),q:ee.List([]),rx5:ee.Image(0)});
 var st=ee.Dictionary(list.iterate(function(obj,acc){
  var img=ee.Image(obj).select('precipitation'); acc=ee.Dictionary(acc);
  var run=ee.Image(acc.get('run')).add(1).where(img.gte(1),0).toInt16();
  var maxrun=ee.Image(acc.get('maxrun')).max(run); var q=ee.List(acc.get('q')).add(img);
  q=ee.List(ee.Algorithms.If(q.length().gt(5),q.slice(q.length().subtract(5)),q));
  var s5=ee.ImageCollection.fromImages(q).sum(), rx5=ee.Image(acc.get('rx5')).max(s5);
  return ee.Dictionary({run:run,maxrun:maxrun,q:q,rx5:rx5});
 },init));
 return rain.addBands(rx1).addBands(ee.Image(st.get('rx5')).rename('rx5day')).addBands(ee.Image(st.get('maxrun')).rename('cdd')).clip(geom).set('year',y);
}
function maskL57(img){
 var qa=img.select('QA_PIXEL'); var m=qa.bitwiseAnd(1<<1).eq(0).and(qa.bitwiseAnd(1<<2).eq(0)).and(qa.bitwiseAnd(1<<3).eq(0)).and(qa.bitwiseAnd(1<<4).eq(0)).and(qa.bitwiseAnd(1<<5).eq(0));
 var sr=img.select(['SR_B3','SR_B4'],['red','nir']).multiply(0.0000275).add(-0.2);
 var st=img.select('ST_B6').multiply(0.00341802).add(149).subtract(273.15).rename('lst');
 return sr.addBands(st).updateMask(m).copyProperties(img,['system:time_start','SPACECRAFT_ID']);
}
function maskL89(img){
 var qa=img.select('QA_PIXEL'); var m=qa.bitwiseAnd(1<<1).eq(0).and(qa.bitwiseAnd(1<<2).eq(0)).and(qa.bitwiseAnd(1<<3).eq(0)).and(qa.bitwiseAnd(1<<4).eq(0)).and(qa.bitwiseAnd(1<<5).eq(0));
 var sr=img.select(['SR_B4','SR_B5'],['red','nir']).multiply(0.0000275).add(-0.2);
 var st=img.select('ST_B10').multiply(0.00341802).add(149).subtract(273.15).rename('lst');
 return sr.addBands(st).updateMask(m).copyProperties(img,['system:time_start','SPACECRAFT_ID']);
}
var landsat=ee.ImageCollection('LANDSAT/LT05/C02/T1_L2').map(maskL57)
 .merge(ee.ImageCollection('LANDSAT/LE07/C02/T1_L2').map(maskL57))
 .merge(ee.ImageCollection('LANDSAT/LC08/C02/T1_L2').map(maskL89))
 .merge(ee.ImageCollection('LANDSAT/LC09/C02/T1_L2').map(maskL89));
function annualLandsat(y){
 y=ee.Number(y); var s=ee.Date.fromYMD(y,1,1),e=s.advance(1,'year'),c=landsat.filterDate(s,e).filterBounds(geom),med=c.median();
 return med.normalizedDifference(['nir','red']).rename('ndvi').addBands(med.select('lst')).clip(geom).set('year',y).set('scene_count',c.size());
}
var validation=ee.FeatureCollection(years.map(function(y){
 var h=annualCHIRPS(y),l=annualLandsat(y);
 var hm=h.reduceRegion({reducer:ee.Reducer.mean(),geometry:geom,scale:5566,maxPixels:1e9,bestEffort:true});
 var lm=l.reduceRegion({reducer:ee.Reducer.mean(),geometry:geom,scale:30,maxPixels:1e10,bestEffort:true});
 return ee.Feature(null,ee.Dictionary({year:y}).combine(hm).combine(lm).set('landsat_scenes',l.get('scene_count')));
}));
print('Compare these basin means with the authoritative annual CSV',validation);
Export.table.toDrive({collection:validation,description:'Tano_Spatial_Raster_Validation',fileFormat:'CSV'});
years.forEach(function(y){
 var h=annualCHIRPS(y),l=annualLandsat(y);
 ['rainfall','rx1day','rx5day','cdd'].forEach(function(b){Export.image.toDrive({image:h.select(b),description:'Tano_'+b+'_'+y,region:geom,scale:5566,maxPixels:1e13,fileFormat:'GeoTIFF'});});
 ['ndvi','lst'].forEach(function(b){Export.image.toDrive({image:l.select(b),description:'Tano_'+b+'_'+y,region:geom,scale:30,maxPixels:1e13,fileFormat:'GeoTIFF'});});
});
