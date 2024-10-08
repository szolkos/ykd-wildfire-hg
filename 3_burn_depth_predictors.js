// ======================================================================================================== //
// Name: 3_burn_depth_predictors
// Author: Scott Zolkos (sgzolkos@gmail.com)
// Created: 2023-01-20
// Revised: 2024-08-20
// Background: Extract remote sensing data for modeling burn depth
// References (see below)
// ======================================================================================================== //

// ========================= //
// ### GLOBAL PARAMETERS ### //
// ========================= //

// COORDINATE SYSTEMS
// EPSG number for output projection (more info: http://spatialreference.org/ref/epsg/); YKD Delta watershed within UTM Zone 3N (https://upload.wikimedia.org/wikipedia/commons/c/cf/UTM_Zones_AK.jpg)
var crs = 'EPSG:32603'; // 'EPSG:32603' = WGS84 / UTM Zone 3N, 'EPSG:4326' = WGS84

// Function to reproject Feature Collection
var fc_rpj = function(feature) {
  var transformed_feature = feature.transform(crs);
  return transformed_feature;
};

// VISUALIZATION
var colorbrewer = require('users/gena/packages:colorbrewer'); // load colorbrewer schemes
var palette_variance = colorbrewer.Palettes.YlOrBr[9];
var l8_vis = {
  bands: ['SR_B4', 'SR_B3', 'SR_B2'],
  min: 0.0,
  max: 0.3,
  gamma: [1.4, 1.4, 1.4],
};
var vis_rgb = {
  min: 0,
  max: 0.2,
  bands: ['SR_B3', 'SR_B2', 'SR_B1']
};
var ndvi_vis = {min: -1, max: 1, palette: ['blue', 'white', 'green']};
var ndwi_vis = {min: -1, max: 1, palette: ['blue', 'white', 'green']};
var nbr_vis = {min: -1, max: 1, palette: ['blue', 'white', 'green', 'yellow', 'orange', 'red', 'black']}; // use if water not masked out
var elevationVis = {
  min: -50,
  max: 80,
  palette: palette_variance,
  //palette: ['0d13d8', '60e1ff', 'ffffff'],
}; // elevation
var slope_vis = {min: 0, max: 7, palette: ['blue', 'white', 'red']}; // slope
var vis_tc = {
  min: [-0.1, 0.0, -0.15],
  max: [0.5, 0.3, 0.0],
  bands: ['brightness', 'greenness', 'wetness']
}; // tasseled-cap (all)
var vis_tc_brightness = {
  min: -0.1, 
  max: 0.5,
  bands: ['brightness'],
  palette: [
    'ffffcc',
    'ffeda0',
    'fed976',
    'feb24c',
    'fd8d3c',
    'fc4e2a',
    'e31a1c',
    'bd0026',
    '800026'
  ]
}; // tasseled-cap (brightness)
var vis_tc_greenness = {
  min: 0.0,
  max: 0.3,
  bands: ['greenness'],
  palette: [
    'ffffe5',
    'f7fcb9',
    'd9f0a3',
    'addd8e',
    '78c679',
    '41ab5d',
    '238443',
    '006837',
    '004529'
  ]
}; // tasseled-cap (greenness)
var vis_tc_wetness = {
  min: -0.15,//-0.1, 
  max: 0.0,//0.1, 
  bands: ['wetness'],
  palette: [
    'fff7fb',
    'ece7f2',
    'd0d1e6',
    'a6bddb',
    '74a9cf',
    '3690c0',
    '0570b0',
    '045a8d',
    '023858'
  ]
}; // tasseled-cap (wetness)

// =================== //
// ### IMPORT DATA ### //
// =================== //

// Sites in YK Delta with burn depth measurements (from Moubarak et al. 2022)
var sites = ee.FeatureCollection("projects/ee-szolkos/assets/NSF_EAR_PF/3_Wildfire_Hg/burn_depth_no_ids");
var sites = sites.select("site", "bd_mean");
// YK Delta watershed boundary from Ludwig et al. 2022 GBC
var ykd_shed = ee.FeatureCollection("projects/ee-szolkos/assets/NSF_EAR_PF/ykd_shed_lg");
// YK Delta surface water polygons from Ludwig et al. 2022 GBC
var water = ee.FeatureCollection("users/luddaludwig/unused_YKD/water_polygons_sr");
var water = water.filterBounds(ykd_shed);
// Fire polygons from Monitoring Trends in Burn Severity (MTBS) & Alaska Interagency Coordination Center (AICC)
var fire_polys_mtbs = ee.FeatureCollection("projects/ee-szolkos/assets/NSF_EAR_PF/3_Wildfire_Hg/mtbs_2015_in_ykd_shed");
var fire_polys = fire_polys_mtbs;
//var fire_polys_aicc = ee.FeatureCollection(fire_polys_aicc).filter('FireYear == "2015"');
//var fire_polys_aicc = ee.FeatureCollection("projects/ee-szolkos/assets/AICC/fire_polys");
// Shapefile of points located wthin Landsat pixel centroids (across YKD watershed), for data extraction
var shed_L8_pts = ee.FeatureCollection("projects/ee-szolkos/assets/NSF_EAR_PF/3_Wildfire_Hg/shed_L8_pts");
var shed_L8_pts_sub = shed_L8_pts.filterBounds(aoi).map(fc_rpj); // Subset ~100 Landsat pixel centroids, for testing prediction
// Shapefile of points located within MTBS 2015 burn scars, for data extraction
//var mtbs_L8_pts = ee.FeatureCollection("projects/ee-szolkos/assets/NSF_EAR_PF/3_Wildfire_Hg/mtbs_L8_pts");
var mtbs_L8_pts = shed_L8_pts.filterBounds(fire_polys);
// ArcticDEMv3 mosaic
var arctic_dem = ee.Image('UMN/PGC/ArcticDEM/V3/2m_mosaic');
// Point data for mean Hg content in the top 30 cm of active layer (organic) soils at field sites that were not located within the perimeters of year 2015 wildfire burn scars
//var ykd_al_hg = ee.FeatureCollection("projects/ee-szolkos/assets/NSF_EAR_PF/3_Wildfire_Hg/ykd_hg_ub_avg");

// ============================ //
// ### LANDSAT BAND METRICS ### //
// ============================ //

// Function to apply required scaling factors (https://www.usgs.gov/faqs/how-do-i-use-scale-factor-landsat-level-2-science-products)
function applyScaleFactors(image) {
  var opticalBands = image.select('SR_B.').multiply(0.0000275).add(-0.2);
  var thermalBands = image.select('ST_B.*').multiply(0.00341802).add(149.0);
  return image.addBands(opticalBands, null, true)
              .addBands(thermalBands, null, true);
}

// IMAGERY & DERIVATIVES
// Landsat 8 Collection 2 Tier 1 Level 2 surface reflectance - pre-2015 wildfire
var l8_prefire = ee.ImageCollection("LANDSAT/LC08/C02/T1_L2")
  .filterDate('2015-06-01','2015-06-21') // Imagery collected 2015-06-16 (according to USGS Earth Explorer)
  .filter(ee.Filter.lt('CLOUD_COVER', 2))
  .filterBounds(ykd_shed);
var l8_postfire = ee.ImageCollection("LANDSAT/LC08/C02/T1_L2")
  .filterDate('2015-08-01','2015-08-31') // NOTE: visual inspection of .filterDate('2016-06-01','2016-08-31') indicates no difference in fire perimeter than '2015-08-01','2015-08-31'
  .filter(ee.Filter.lt('CLOUD_COVER', 2))
  .filterBounds(ykd_shed);
// Note: Moubarak et al. 2022 defined pre- and post-fire date ranges as 6/01–8/31 (pre = 2014, post = 2016), but there appears to be insufficient cloud-free imagery. Hence, dates above collect imagery immediately before & after fire in 2015.

// Apply scaling factors
var l8_prefire = l8_prefire.map(applyScaleFactors).mean().clip(ykd_shed);
var l8_postfire = l8_postfire.map(applyScaleFactors).mean().clip(ykd_shed);

// Subset desired bands & define projection
var l8_prefire = l8_prefire.select('SR_B2','SR_B3','SR_B4','SR_B5','SR_B6','SR_B7');
var l8_prefire = l8_prefire.reproject(crs, null, 30); // reproject from 2m to 30m and define projected coordinate system
var l8_postfire = l8_postfire.select('SR_B2','SR_B3','SR_B4','SR_B5','SR_B6','SR_B7');
var l8_postfire = l8_postfire.reproject(crs, null, 30);
//var l8_b234 = l8_postfire.select('SR_B4','SR_B3','SR_B2');
//var l8_b234 = l8_b234.reproject(crs, null, 30);
var l8_b2 = l8_prefire.select('SR_B2');
var l8_b2 = l8_b2.reproject(crs, null, 30);
var l8_b3 = l8_prefire.select('SR_B3');
var l8_b3 = l8_b3.reproject(crs, null, 30);
var l8_b4 = l8_prefire.select('SR_B4');
var l8_b4 = l8_b4.reproject(crs, null, 30);

// Calculate pre- & post-fire NBR = (NIR-SWIR2)/(NIR+SWIR2)
// https://developers.google.com/earth-engine/apidocs/ee-image-normalizeddifference
var nbr_pre = l8_prefire.normalizedDifference(['SR_B5', 'SR_B7']).rename('preNBR'); // B5 (0.85-0.88 µm) = NIR, B7 (2.11-2.29 µm) = SWIR2
var nbr_post = l8_postfire.normalizedDifference(['SR_B5', 'SR_B7']).rename('postNBR');

// dNBR = NBRprefire - NBRpostfire * 1000
var dnbr = nbr_pre.subtract(nbr_post).rename('dNBR');//.multiply(1000);
//var dnbr = dnbr.updateMask(dnbr.gte(0));

// RBR = dNBR / (NBRprefire + 1.001) (Parks et al. 2014 RS)
var rbr = dnbr.divide(nbr_pre.add(1.001)).rename('RBR');

// RdnBR = dNBR / (NBRprefire + 1.001) (Miller and Thode 2007 RSE)
var rdnbr = dnbr.divide(((nbr_pre.divide(1000)).abs()).sqrt()).rename('RdNBR');

// NDVI (NIR-Red)/(NIR+Red) (b2,3,4,5 = b,g,r,nir)
var ndvi = l8_prefire.normalizedDifference(['SR_B5', 'SR_B4']).rename('NDVI');

// NDWI (G-NIR)/(G+NIR) | https://developers.google.com/earth-engine/tutorials/tutorial_api_06
var ndwi = l8_prefire.normalizedDifference(['SR_B3', 'SR_B5']).rename('NDWI'); // B5 (0.85-0.88 µm) = NIR

// MNDWI (https://www.space4water.org/taxonomy/term/1246)

// NDII1 (NIR-SWIR1)/(NIR+SWIR1) (Hardisky et al. 1983 PERS)
var ndii1 = l8_prefire.normalizedDifference(['SR_B5', 'SR_B6']).rename('NDII1'); // B5 (0.85-0.88 µm) = NIR, B6 (1.57-1.65 µm) = SWIR1

// NDII2 (NIR-SWIR2)/(NIR+SWIR2) 
var ndii2 = l8_prefire.normalizedDifference(['SR_B5', 'SR_B7']).rename('NDII2'); // B5 (0.85-0.88 µm) = NIR, B7 (2.11-2.29 µm) = SWIR2


// ============================ //
// ### TASSELED-CAP INDICES ### //
// ============================ //

// Interpretation of RGB visualization of Brightness, Greenness, and Wetness:
// https://desktop.arcgis.com/en/arcmap/10.3/manage-data/raster-and-images/tasseled-cap-transformation.htm

// FUNCTIONS
/**
 * Description of function
 * @param   {param_type}   param_name  param_description
 * @return  {return_type}              return_description
 */

/**
 * Creates a Tasseled Cap image for Landsat 5
 * @param   {ee.Image}  image  Landsat 5 image
 * @return  {ee.Image}         Tasseled Cap image with 6 bands: 
 *                               'brightness', 'greenness', 'wetness',
 *                               'fourth', 'fifth', sixth'
 */ 
function tasseled_cap_L5(image) {
  // Define array of Tasseled Cap coefficients
  var coefficients = ee.Array([
    [  0.3037,  0.2793,  0.4743,  0.5585,  0.5082,  0.1863 ],
    [ -0.2848, -0.2435, -0.5436,  0.7243,  0.0840, -0.1800 ],
    [  0.1509,  0.1973,  0.3279,  0.3406, -0.7112, -0.4572 ],
    [ -0.8242,  0.0849,  0.4392, -0.0580,  0.2012, -0.2768 ],
    [ -0.3280,  0.0549,  0.1075,  0.1855, -0.4357,  0.8085 ],
    [  0.1084, -0.9022,  0.4120,  0.0573, -0.0251,  0.0238 ]
  ]);
  
  // Select bands for use in Tasseled Cap
  var image_bands_tc = image.select(['B1', 'B2', 'B3', 'B4', 'B5', 'B7']);
  
  // Create 1-D array image (vector of length 6 for all bands per pixel)
  var array_image_1d = image_bands_tc.toArray();

  // Create 2-D array image (6x1 matrix for all bands per pixel) from 1-D array
  var array_image_2d = array_image_1d.toArray(1);
  
  // Get a multi-band image with TC-named bands 
  // Matrix multiplication: 6x6 times 6x1
  var components_image = ee.Image(coefficients)
    .matrixMultiply(array_image_2d)
    // Remove extra dimensions
    .arrayProject([0])
    // Convert to regular image
    .arrayFlatten([['brightness', 'greenness', 'wetness', 'fourth', 'fifth', 'sixth']]);
    
  return components_image;
}

// Tasseled-Cap transformation for Landsat 8 OLI surface reflectance data
// Coefficients obtained from Zhai et al. 2022 RSE
function tasseled_cap_L8(image) {
  // Define array of Tasseled Cap coefficients
  var coefficients = ee.Array([
    [  0.3690,  0.4271,  0.4689,  0.5073,  0.3824,  0.2406 ], // Brightness
    [ -0.2870, -0.2685, -0.4087,  0.8145,  0.0637, -0.1052 ], // Greenness
    [  0.0382,  0.2137,  0.3536,  0.2270, -0.6108, -0.6351 ] // Wetness
  ]);
  
  // Select bands for use in Tasseled Cap: B2 (blue), B3 (green), B4 (red), B5 (NIR), B6 (SWIR1), B7 (SWIR2)
  var image_bands_tc = image.select(['SR_B2', 'SR_B3', 'SR_B4', 'SR_B5', 'SR_B6', 'SR_B7']);
  
  // Create 1-D array image (vector of length 3 for all bands per pixel)
  var array_image_1d = image_bands_tc.toArray();

  // Create 2-D array image (3x1 matrix for all bands per pixel) from 1-D array
  var array_image_2d = array_image_1d.toArray(1);
  
  // Get a multi-band image with TC-named bands 
  // Matrix multiplication: 6x3 times 6x1
  var components_image = ee.Image(coefficients)
    .matrixMultiply(array_image_2d)
    // Remove extra dimensions
    .arrayProject([0])
    // Convert to regular image
    .arrayFlatten([['brightness', 'greenness', 'wetness']]);
    
  return components_image;
}

// Create pre-fire tasseled-cap indices
var ykd_tc = tasseled_cap_L8(l8_prefire);
var ykd_tc = ykd_tc.reproject(crs, null, 30);


// ========================= //
// ### ELEVATION METRICS ### //
// ========================= //

// Elevation (ArcticDEM v3 2m mosaic, digital surface model)
var elev = arctic_dem.select('elevation');
var elev = elev.reproject(crs, null, 30).clip(ykd_shed); // reproject from 2m to 30m and define projected coordinate system

// Slope (units = degrees, range is [0,90])
var slp = ee.Terrain.slope(elev);

// Aspect (units = degrees, where 0=N, 90=E, 180=S, 270=W)
var asp = ee.Terrain.aspect(elev);

// Rugosity (measure of small-scale variations of amplitude in surface height)
var rugos = arctic_dem.reduceNeighborhood({
  reducer: ee.Reducer.stdDev(), //ee.Reducer.variance
  kernel: ee.Kernel.square(15),
});
var rugos = rugos.reproject(crs, null, 30).clip(ykd_shed); 


// =================================== //
// ### PREPARE DATA FOR EXTRACTION ### //
// =================================== //

// Concatenate all data above as multiple bands of an image (merge predictor variables into one ImageCollection)
var images = ee.Image.cat(
  l8_b2,
  l8_b3,
  l8_b4,
  dnbr,
  //rbr,
  //rdnbr,
  ndvi,
  ndwi,
  //ndii1,
  //ndii2,
  ykd_tc.select('brightness'),
  ykd_tc.select('greenness'), 
  ykd_tc.select('wetness'),
  elev,
  slp,
  rugos)
  //asp)
  .set('system:time_start', ee.Date('2000-01-01').millis()); // Computed images do not have a 'system:time_start' property; add one based on when the data were collected.

// Wrap the single image in an ImageCollection for use in the zonalStats function.
var ic = ee.ImageCollection([images]);
//var ic = ic.clip(fire_polys);

// Mask out water
var mask = ee.Image.constant(1).clip(water.geometry()).mask().not();
var preds = ic.map(function(image){return image.updateMask(mask)});

// Define parameters to extract data at BD mmnt sites using the below functions ('zonalStats') and for Landsat pixel centroids ('fill')
var params = {
  bands: [0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16],
  bandsRename: ['l8_b2', 
                'l8_b3',
                'l8_b4',
                'dnbr',
                //'rbr',
                //'rdnbr',
                'ndvi',
                'ndwi',
                //'ndii1',
                //'ndii2',
                'tc_bright',
                'tc_green',
                'tc_wet',
                'elev',
                'slope',
                'rugosity']//,
                //'asp']
};

// ========================================= //
// ### BURN DEPTH (BD) MODELING STEP 1:  ###
// ### EXTRACT DATA AT BD MMNT LOCATIONS ### //
// ### FOR MODEL DEVELOPMENT (PREDICTORS)###
// ========================================= //
// Adapted from: https://developers.google.com/earth-engine/tutorials/community/extract-raster-values-for-points
// Buffer points
function bufferPoints(radius, bounds) {
  return function(pt) {
    pt = ee.Feature(pt);
    return bounds ? pt.buffer(radius).bounds() : pt.buffer(radius);
  };
}

// Zonal stats (see: https://developers.google.com/earth-engine/tutorials/community/extract-raster-values-for-points)
function zonalStats(ic, fc, params) {
  // Initialize internal params dictionary.
  var _params = {
    reducer: ee.Reducer.mean(),
    scale: null,
    crs: null,
    bands: null,
    bandsRename: null,
    imgProps: null,
    imgPropsRename: null,
    datetimeName: 'datetime',
    datetimeFormat: 'YYYY-MM-dd HH:mm:ss'
  };

  // Replace initialized params with provided params.
  if (params) {
    for (var param in params) {
      _params[param] = params[param] || _params[param];
    }
  }

  // Set default parameters based on an image representative.
  var imgRep = ic.first();
  var nonSystemImgProps = ee.Feature(null)
    .copyProperties(imgRep).propertyNames();
  if (!_params.bands) _params.bands = imgRep.bandNames();
  if (!_params.bandsRename) _params.bandsRename = _params.bands;
  if (!_params.imgProps) _params.imgProps = nonSystemImgProps;
  if (!_params.imgPropsRename) _params.imgPropsRename = _params.imgProps;

  // Map the reduceRegions function over the image collection.
  var results = ic.map(function(img) {
    // Select bands (optionally rename), set a datetime & timestamp property.
    img = ee.Image(img.select(_params.bands, _params.bandsRename))
      .set(_params.datetimeName, img.date().format(_params.datetimeFormat))
      .set('timestamp', img.get('system:time_start'));

    // Define final image property dictionary to set in output features.
    var propsFrom = ee.List(_params.imgProps)
      .cat(ee.List([_params.datetimeName, 'timestamp']));
    var propsTo = ee.List(_params.imgPropsRename)
      .cat(ee.List([_params.datetimeName, 'timestamp']));
    var imgProps = img.toDictionary(propsFrom).rename(propsFrom, propsTo);

    // Subset points that intersect the given image.
    var fcSub = fc.filterBounds(img.geometry());

    // Reduce the image by regions.
    return img.reduceRegions({
      collection: fcSub,
      reducer: _params.reducer,
      scale: _params.scale,
      crs: _params.crs
    })
    // Add metadata to each feature.
    .map(function(f) {
      return f.set(imgProps);
    });
  }).flatten().filter(ee.Filter.notNull(_params.bandsRename));

  return results;
}

// Buffer burn depth mmnt site points
var sites_bfr = sites.map(bufferPoints(2, false));

// Buffer Landsat 8 pixel centroids within MTBS 2015 burn scars, to extract values for predictive modeling of burn depth in R
//var sites_bfr = mtbs_L8_pts.map(bufferPoints(2, false));

// Extract zonal statistics per point per image.
//var bd_site_stats = zonalStats(ic, sites_bfr, params);


// ========================================= //
// ### BURN DEPTH (BD) MODELING STEP 2:  ###
// ### EXTRACT DATA AT LANDSAT PIXEL     ### //
// ### CENTROIDS FOR PREDICTING BD       ###
// ========================================= //
// # Adapted from: https://stackoverflow.com/questions/42742742/extract-pixel-values-by-points-and-convert-to-a-table-in-google-earth-engine

// GET PIXEL CENTROID COORDINATES
var latlon = ee.Image.pixelLonLat().reproject(crs, null, 30);

// Set region within which pixel centroids will be determined
var bounding_geom = fire_polys;

// Put coordinates in a list
var coords = latlon.select(['longitude', 'latitude']).reduceRegion({
  reducer: ee.Reducer.toList(),
  geometry: bounding_geom,
  scale: 30,
  maxPixels: 10e7
});

// Get lat & lon
var lat = ee.List(coords.get('latitude'));
var lon = ee.List(coords.get('longitude'));

// Zip lat & lon 
// (e.g.: zip([1, 3],[2, 4]) --> [[1, 2], [3,4]])
var point_list = lon.zip(lat);

// Create points from coordinates
var centroids = ee.Geometry.MultiPoint(point_list);

// Convert centroids from Multipoint to FeatureCollection
var centroids = ee.FeatureCollection(centroids.coordinates().map(function(p){
  var point = ee.Feature(ee.Geometry.Point(p), {});
  return point}));

// EXTRACT VALUES FROM IMAGE COLLECTION AT EACH POINT
// Empty Collection to fill
var ft = ee.FeatureCollection(ee.List([]));

var fill = function(img, ini) {
  // type cast
  var inift = ee.FeatureCollection(ini);

  // gets the values for the points in the current img
  var ft2 = img.reduceRegions(centroids, ee.Reducer.first(), 30);

  // gets the date of the img
  //var date = img.date().format()

  // writes the date in each feature
  //var ft3 = ft2.map(function(f){return f.set("date", date)})

  // merges the FeatureCollections
  return inift.merge(ft2)};

// Iterates over the ImageCollection
var pixel_centroid_stats = ee.FeatureCollection(ic.iterate(fill, ft));


// ====================== //
// ### DATA SUMMARIES ### //
// ====================== //

print("Burn depth mmnt points (transect average):", sites);
//print("Burn depth mmnt points (transect average) (+0.1 m buffer):", sites_bfr);
//print("Predictor variables:", preds);
//print("YKD burn depth mmnt site stats: ", pixel_centroid_stats);
//print("Centroids (correct, from GEE) projections:", centroids.filterBounds(aoi));
//print("Centroids (off, from ArcGIS Pro) projections:", shed_L8_pts_sub);
//print("Burn depth (ranger random forest):", burn_depth_ranger);

// =================== //
// ### EXPORT DATA ### //
// =================== //

/*
Export.table.toDrive({
  collection: pixel_centroid_stats,
  //description: 'burn_depth_preds_gee_2023_01_22', // name of task and by default, name of csv file (no spaces allowed)
    description: 'soc_preds_from_gee_2023_10_06',
  folder: 'Geospatial', // replace with name of folder in your Google Drive
  fileFormat: 'CSV'
});
*/
