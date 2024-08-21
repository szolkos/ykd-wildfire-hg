// ======================================================================================================== //
// Name: 4_soc_predictors
// Author: Scott Zolkos (sgzolkos@gmail.com)
// Created: 2023-01-05
// Revised: 2024-08-20
// Background: Extract remote sensing data for modeling SOC stores
// References (see below)
// ======================================================================================================== //

// ========================= //
// ### GLOBAL PARAMETERS ### //
// ========================= //

// IMPORT DATA
var ykd_shed = ee.FeatureCollection("projects/ee-szolkos/assets/NSF_EAR_PF/ykd_shed_lg");
var water = ee.FeatureCollection("users/luddaludwig/unused_YKD/water_polygons_sr");
var water = water.filterBounds(ykd_shed);
var fire_polys_mtbs = ee.FeatureCollection("projects/ee-szolkos/assets/MTBS/mtbs_2015_in_ykd_shed");
var fire_polys = fire_polys_mtbs;
var mtbs_L8_pts = ee.FeatureCollection("projects/ee-szolkos/assets/NSF_EAR_PF/3_Wildfire_Hg/mtbs_L8_pts");
var shed_L8_pts = ee.FeatureCollection("projects/ee-szolkos/assets/NSF_EAR_PF/3_Wildfire_Hg/shed_L8_pts");
//var ykd_al_hg = ee.FeatureCollection("projects/ee-szolkos/assets/NSF_EAR_PF/3_Wildfire_Hg/ykd_hg_ub_avg");
var soc_site = ee.FeatureCollection("projects/ee-szolkos/assets/NSF_EAR_PF/3_Wildfire_Hg/soc_field_sites_2023_10_06");
var hg_site = ee.FeatureCollection("projects/ee-szolkos/assets/NSF_EAR_PF/3_Wildfire_Hg/hg_field_sites_2023_10_06");

// COORDINATE SYSTEMS
// EPSG number for output projection (more info: http://spatialreference.org/ref/epsg/); YKD Delta watershed within UTM Zone 3N (https://upload.wikimedia.org/wikipedia/commons/c/cf/UTM_Zones_AK.jpg)
var crs = 'EPSG:32603'; // 'EPSG:32603' = WGS84 / UTM Zone 3N, 'EPSG:4326' = WGS84

// VISUALIZATION
// Load colorbrewer schemes
var colorbrewer = require('users/gena/packages:colorbrewer'); // load colorbrewer schemes
var palette_variance = colorbrewer.Palettes.YlOrBr[9];
var l8_vis = {
  bands: ['SR_B4', 'SR_B3', 'SR_B2'],
  min: 0.0,
  max: 0.3,
  gamma: [1.4, 1.4, 1.4],
}; // Define Landsat 8 vis parameters
var vis_rgb = {
  min: 0,
  max: 0.2,
  bands: ['SR_B3', 'SR_B2', 'SR_B1']
};
var ndvi_vis = {min: -1, max: 1, palette: ['blue', 'white', 'green']};
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

// Landsat 8 Collection 2 Tier 1 Level 2 surface reflectance - pre-2015 wildfire
var l8 = ee.ImageCollection("LANDSAT/LC08/C02/T1_L2")
  .filterDate('2015-06-01','2015-06-21') // pre-fire
  .filter(ee.Filter.lt('CLOUD_COVER', 2))
  .filterBounds(ykd_shed);

// Apply scaling factors
var l8 = l8.map(applyScaleFactors).mean().clip(ykd_shed);

// Subset desired bands
var l8_b234 = l8.select('SR_B4','SR_B3','SR_B2');
var l8_b234 = l8_b234.reproject(crs, null, 30); // reproject from 2m to 30m and define projected coordinate system
var l8_b2 = l8_b234.select('SR_B2');
var l8_b3 = l8_b234.select('SR_B3');
var l8_b4 = l8_b234.select('SR_B4');

// NDVI (NIR-Red)/(NIR+Red) (b2,3,4,5 = b,g,r,nir)
var ndvi = l8.normalizedDifference(['SR_B5', 'SR_B4']).rename('NDVI');

// NDWI (G-NIR)/(G+NIR) | https://developers.google.com/earth-engine/tutorials/tutorial_api_06
var ndwi = l8.normalizedDifference(['SR_B3', 'SR_B5']).rename('NDWI');


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

// Apply required scaling factors to Landsat SR imagery (https://www.usgs.gov/faqs/how-do-i-use-scale-factor-landsat-level-2-science-products)
function applyScaleFactors(image) {
  var opticalBands = image.select('SR_B.').multiply(0.0000275).add(-0.2);
  var thermalBands = image.select('ST_B.*').multiply(0.00341802).add(149.0);
  return image.addBands(opticalBands, null, true)
              .addBands(thermalBands, null, true);
}

// Create tasseled-cap indices
var ykd_tc = tasseled_cap_L8(l8);
var ykd_tc = ykd_tc.reproject(crs, null, 30);


// ========================= //
// ### ELEVATION METRICS ### //
// ========================= //

// Elevation (ArcticDEM v3 2m mosaic, digital surface model)
var arctic_dem = ee.Image('UMN/PGC/ArcticDEM/V3/2m_mosaic').clip(ykd_shed);
var elev = arctic_dem.select('elevation').clip(ykd_shed);
//var elev = elev.reproject(crs, null, 30); // reproject from 2m to 30m and define projected coordinate system
// Resample from 2 m to 30 m using bilinear interpolation and define projected coordinate system
var elev = elev
  .resample('bilinear')
  .reproject(crs, null, 30);
  
// Slope (units = degrees, range is [0,90])
var slp = ee.Terrain.slope(elev);

// Aspect (units = degrees, where 0=N, 90=E, 180=S, 270=W)
var asp = ee.Terrain.aspect(elev);

// Rugosity (measure of small-scale variations of amplitude in surface height)
var rugos = arctic_dem.reduceNeighborhood({
  reducer: ee.Reducer.stdDev(), //ee.Reducer.variance
  kernel: ee.Kernel.square(15),
});
print("Rugosity: ", rugos);


// =================================== //
// ### PREPARE DATA FOR EXTRACTION ### //
// =================================== //

// Concatenate all data above as multiple bands of an image (merge predictor variables into one ImageCollection)
var images = ee.Image.cat(
  l8_b2,
  l8_b3,
  l8_b4,
  ndvi,
  ndwi,
  ykd_tc.select('brightness'),
  ykd_tc.select('greenness'), 
  ykd_tc.select('wetness'),
  slp,
  asp,
  elev,
  rugos)
  .set('system:time_start', ee.Date('2000-01-01').millis()); // Computed images do not have a 'system:time_start' property; add one based on when the data were collected.

// Wrap the single image in an ImageCollection for use in the zonalStats function.
var ic = ee.ImageCollection([images]);


// ======================================= //
// ### EXTRACT DATA AT SAMPLING POINTS ### //
// ======================================= //
// https://developers.google.com/earth-engine/tutorials/community/extract-raster-values-for-points
// Select to extract data at sites with data for SOC (n=51) or Hg (n=50)
var site_data = hg_site; // soc_site or hg_site

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

// Buffer sampling points
var site_data_bfr = site_data.map(bufferPoints(2, false));

// Define parameters for the zonalStats function.
var params = {
  bands: [0,1,2,3,4,5,6,7,8,9,10,11],
  bandsRename: ['l8_b2', 'l8_b3', 'l8_b4', 'ndvi', 'ndwi', 'tc_bright', 'tc_green', 'tc_wet', 'slope', 'asp', 'elev', 'rugosity']
};

// Extract zonal statistics per point per image.
var pts_stats = zonalStats(ic, site_data_bfr, params);
print("YKD Hg sampling site stats: ", pts_stats);


// ====================== //
// ### DATA SUMMARIES ### //
// ====================== //

print("Hg sampling points:", site_data);
print("Hg sampling points (+0.1 m buffer):", site_data_bfr);
print('YK Delta Tasseled Cap:', ykd_tc);
//print('YK Delta pre-fire mean:', l8);
//print("MTBS L8 points: ", mtbs_L8_pts.first());


// =================== //
// ### EXPORT DATA ### //
// =================== //

/*
Export.table.toDrive({
  collection: pts_stats,
  description: 'preds_for_soil_hg_modeling_2023_10_06', // name of task and by default, name of csv file (no spaces allowed)
  folder: 'Geospatial', // replace with name of folder in your Google Drive
  fileFormat: 'CSV'
});
*/
