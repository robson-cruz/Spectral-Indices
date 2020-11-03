/*=====================================================================================
          
           INDICES ESPECTRAIS PARA REALÇAR VEGETAÇÃO NO ESTADO DO PARÁ
          
  =====================================================================================
*/

/*
* NDVI (Normalized Difference Vegetation Index)
* EVI (Enhanced Vegetation Index)
* SAVI (Soil Adjusted Vegetation Index)
* NDFI (Normalized Difference Fraction Index)
* DETEX (INPE)
*/


// Load State of Pará.
//var pa = ee.FeatureCollection('users/nsrditec/IBGE/municipios_pa');

// Load Municipal Mesh of Pará.
var municipio = ee.FeatureCollection('users/nsrditec/IBGE/municipios_pa')
  // Select Municipality
  .filterMetadata('NM_MUNICIP', 'equals', 'TERRA SANTA');

// Load Forest Manage Unit (UMF) - Concession Forest.
var UMF = ee.FeatureCollection('users/nsrditec/IBAMA/umf_concessao_florestal_pa')
  .filterMetadata('flona', 'equals', 'Altamira')
  .filterMetadata('umf', 'equals', 'UMF-III');  

print('UMF:', UMF);

var umfGeom = UMF.geometry();
var bufferUmf = umfGeom.buffer(5000);

// Load the polygon of the AOI.
var UPA = ee.FeatureCollection('users/nsrditec/IBAMA/UPA_PA')
  .filterMetadata('Ano_POA', 'equals', '2020')
  .filterMetadata('UMF', 'equals', '3')
  .filterMetadata('FLONA', 'equals', 'Altamira');

print(UPA);

var upaGeom = UPA.geometry();
var bufferUpa = upaGeom.buffer(5000);

/*
**********************************************************************************
*     Functions for Fusion Image Collections of LANDSAT Satellites, TM, ETM+     *
*     and OLI sensors.                                                           *
**********************************************************************************
*/

// Define coefficients supplied by Roy et al. (2016) for translating ETM+
// surface reflectance to OLI surface reflectance.
var coefficients = {
  itcps: ee.Image.constant([0.0003, 0.0088, 0.0061, 0.0412, 0.0254, 0.0172])
        .multiply(10000),
  slopes: ee.Image.constant([0.8474, 0.8483, 0.9047, 0.8462, 0.8937, 0.9071])
};

// Define function to get and rename bands of interest from OLI.
function renameOLI(img){
  return img.select(
		['B2', 'B3', 'B4', 'B5', 'B6', 'B7', 'pixel_qa'],
		['Blue', 'Green', 'Red', 'NIR', 'SWIR1', 'SWIR2', 'pixel_qa']);
}

// Define function to get and rename bands of interest from TM.
function renameTM(img){
  return img.select(
		['B1', 'B2', 'B3', 'B4', 'B5', 'B7', 'pixel_qa'],
		['Blue', 'Green', 'Red', 'NIR', 'SWIR1', 'SWIR2', 'pixel_qa']);
}

// Define function to apply harmonization transformation.
function tm2oli(img){
  return img.select(['Blue', 'Green', 'Red', 'NIR', 'SWIR1', 'SWIR2'])
    .multiply(coefficients.slopes)
    .add(coefficients.itcps)
    .round()
    .toShort()
    .addBands(img.select('pixel_qa'));
}

// Define function to mask out clouds and cloud shadows.
function fmask(img){
  var cloudShadowBitMask = 1 << 3;
  var cloudsBitMask = 1 << 5;
  var qa = img.select('pixel_qa');
  var mask = qa.bitwiseAnd(cloudShadowBitMask).eq(0)
            .and(qa.bitwiseAnd(cloudsBitMask).eq(0));
  return img.updateMask(mask);
}

// Define Function NDVI - Normalized Difference Vegetation Index.
var NDVI = function(img){
  return img.expression(
    '(NIR - Red) / (NIR + Red)',
    {
      'NIR': img.select('NIR'),
      'Red': img.select('Red')
    }).rename('NDVI').set('system:time_start', img.get('system:time_start'));
};

// Define Function EVI - Enhanced Vegetation Index.
var EVI = function(img){
  return img.expression(
    'G * ((NIR - Red) / (L + NIR + C1 * Red - C2 * Blue))', 
    {
      'G': 2.5,
      'L': 1.0,
      'C1': 6.0,
      'C2': 7.5,
      'NIR': img.select('NIR').divide(10000),
      'Red': img.select('Red').divide(10000),
      'Blue': img.select('Blue').divide(10000)
    }).rename('EVI').set('system:time_start', img.get('system:time_start'));
  };

// Define Function NDWI - Normalized difference Water Index.
var NDWI = function(img){
  return img
    .normalizedDifference(['Green', 'NIR'])
    .rename('NDWI')
    .set('system:time_start', img.get('system:time_start'));
};

// Define Function SAVI - Soil Adjusted Vegetation Index.
var SAVI = function(img){
  return img.expression(
      '(1 + L) * (NIR - Red) / (NIR + Red + L)',
      {
        'NIR': img.select('NIR').divide(10000),
        'Red': img.select('Red').divide(10000),
        'L': 0.5
      })
      .rename('SAVI').set('system:time_start', img.get('system:time_start'));
};

// Define Function NDFI - Normalized Difference Fraction Index (SOUZA Jr et al., 2005)
// Intact Forest: NDFI > 0.75
// Canopy Damage: 0 <= NDFI <= 0.75

function NDFI(img){ // Souza  Jr.  et  al  (2005)
  var GV = [500, 900, 400, 6100, 3000, 1000];
  var Shade = [0, 0, 0, 0, 0, 0];
  var NPV = [1400, 1700, 2200, 3000, 5500, 3000];
  var Soil = [2000, 3000, 3400, 5800, 6000, 5800];
  var Cloud = [9000, 9600, 8000, 7800, 7200, 6500];
  var unmixed = img
    .select(['Blue', 'Green', 'Red', 'NIR', 'SWIR1', 'SWIR2'])
    .unmix({
      endmembers: [GV, Shade, NPV, Soil, Cloud],
      sumToOne: true,
      nonNegative: true
    })
    .rename(['GV', 'Shade', 'NPV', 'Soil', 'Cloud']);
  return unmixed
    .expression(
      '((i.GV / (1 - i.Shade)) - (i.NPV + i.Soil)) / ((i.GV / (1 - i.Shade)) + i.NPV + i.Soil)',
      {
        'i': unmixed
      }
    ) 
    .rename('NDFI').set('system:time_start', img.get('system:time_start'));
}


// Define Function DETEX
// DETEX = GV * soil /  
// Separação das frações imagem mais recente
var DETEX = function (img){
  // Define Endmembers
  var soil = [1875, 2263, 2763]; // Soil (Solo).
  var GV = [221, 3480.5, 1559.5]; // GV (Vegetação Verde).
  var shade = [0, 0, 0]; // Shade (Sombra/Água).
  // Uminxing Data
  var unmixed = img
    .select(['Red','NIR','SWIR1'])
    .unmix({endmembers:[soil, GV, shade]})
    .multiply(100)
    .max(0)
    .byte()
    .rename(['soil', 'GV', 'shade']);
    // Cálculo do Detex
  return unmixed
    .expression(
      '(120 * (i.soil / i.GV) + 0)',
      {
        'i': unmixed,
        'soil': unmixed.select('soil'),
        'GV': unmixed.select('GV')
      })
    .rename('DETEX').set('system:time_start', img.get('system:time_start'));
};

// Define function to prepare OLI images.
function prepOLI(img){
  var orig = img;
  img = renameOLI(img);
  img = fmask(img);
  return ee.Image(img.copyProperties(orig, orig.propertyNames()));
}
  
// Define function to prepare TM images.
function prepTM(img){
  var orig = img;
  img = renameTM(img);
  img = fmask(img);
  img = tm2oli(img);
  return ee.Image(img.copyProperties(orig, orig.propertyNames()));
}

// Define function to resample img by bicubic.
function resampleImg(img){
  return img.resample('bicubic'); // bicubic, bilinear, nearest-neighbor
}

// Get Landsat surface reflectance collections for OLI, ETM+ and TM sensors
var oliCol = ee.ImageCollection('LANDSAT/LC08/C01/T1_SR');
var etmCol= ee.ImageCollection('LANDSAT/LE07/C01/T1_SR');
var tmCol= ee.ImageCollection('LANDSAT/LT05/C01/T1_SR');

// Define a first collection filter.
var colFilter_1 = ee.Filter.and(
  ee.Filter.bounds(UMF),
  ee.Filter.date('2019-01-01', '2019-12-31'),
  ee.Filter.lt('CLOUD_COVER', 16),
  ee.Filter.lt('GEOMETRIC_RMSE_MODEL', 10),
  ee.Filter.or(
    ee.Filter.eq('IMAGE_QUALITY', 9),
    ee.Filter.eq('IMAGE_QUALITY_OLI', 9)));

// Define a second collection filter.
var colFilter_2 = ee.Filter.and(
  ee.Filter.bounds(UMF),
  ee.Filter.date('2020-01-01', '2020-10-20'),
  ee.Filter.lt('CLOUD_COVER', 16),
  ee.Filter.lt('GEOMETRIC_RMSE_MODEL', 10),
  ee.Filter.or(
    ee.Filter.eq('IMAGE_QUALITY', 9),
    ee.Filter.eq('IMAGE_QUALITY_OLI', 9)));

// Filter collections and prepare them for merging - first collection.
var oliCol_1 = oliCol.filter(colFilter_1).map(prepOLI).map(resampleImg);
var etmCol_1 = etmCol.filter(colFilter_1).map(prepTM).map(resampleImg);
var tmCol_1 = tmCol.filter(colFilter_1).map(prepTM).map(resampleImg);

// Filter collections and prepare them for merging - second collection.
var oliCol_2 = oliCol.filter(colFilter_2).map(prepOLI);
var etmCol_2 = etmCol.filter(colFilter_2).map(prepTM);
var tmCol_2 = tmCol.filter(colFilter_2).map(prepTM);

// Merge the collections.
var col_1 = oliCol_1
  .merge(etmCol_1)
  .merge(tmCol_1);

var col_2 = oliCol_2
  .merge(etmCol_2)
  .merge(tmCol_2);

// Get collectiona information 
var infColMedian1 = col_1.sort('system:time_start', false);
print('SÉRIE TEMPORAL MEDIANA Col_1:', infColMedian1);
var infColMedian2 = col_2.sort('system:time_start', false);
print('SÉRIE TEMPORAL MEDIANA Col_2:', infColMedian2);

// Reduce the ImageCollection to intra-annual median.
function medianAnual(img){
  return ee.ImageCollection(img
    .reduce(ee.Reducer.median()))
    .mosaic();
} 

function INDEXES(img){
  return img
  .map(EVI)
  .reduce(ee.Reducer.median())
  .clipToCollection(municipio);
}

var test = INDEXES(col_1);
var medianAnual_1 = medianAnual(col_1).clip(UMF);
var medianAnual_2 = medianAnual(col_2).clip(UMF);

/*
*********************************************************************************
*                 RGB Display Parameters and Spectral Indexes                   *
*********************************************************************************
*/
var visNDVI = {
  min: 0.2,
  max: 0.9,
  gamma: 1.4
  /*palette: [
      'FFFFFF', 'CE7E45', 'DF923D', 'F1B555', 'FCD163', '99B718',
      '74A901', '66A000', '529400', '3E8601', '207401', '056201',
      '004C00', '023B01', '012E01', '011D01', '011301'
  ]*/
};

var visNDWI = {
  min: 0,
  max: 1,
  palette: ['#D3D3D3', '#00FFFF'],
};

var visNDFI = {
  min: 0.4,
  max: 1,
  gamma: 1
};

var rgbMedian = {
  bands: ['SWIR1_median', 'NIR_median', 'Red_median'],
  min: 0,
  max: 10000,
  gamma: 1.4
};

// Export Index or Mosaic Images to Google Drive
Export.image.toDrive({
  image: col_2.map(NDFI).median(),
  description: 'NDFI_8OLI_2020_UMF_1B_Saraca_taquera',
  folder: 'img',
  region: bufferUmf,
  scale: 30,
  maxPixels: 10e9,
  fileFormat: 'GeoTiff'
});

// Export Index Images Metadata to Google Drive
Export.table.toDrive({
  collection: col_2, 
  description: 'NDFI_images_metadata_2020_umf_1B_Saraca_taquera', 
  folder: 'img', 
  fileFormat: 'CSV'
  }); 
  
/*
*********************************************************************************
*--------------------------------SHOW RESULTS-----------------------------------*
*********************************************************************************
*/
Map.centerObject(bufferUmf, 12);
//Map.addLayer(col_1.map(EVI).median().clip(municipio),visNDVI,'EVI');
//Map.addLayer(col_2.map(NDVI).median().clip(municipio),visNDVI,'NDVI');
Map.addLayer(medianAnual_1, rgbMedian,'RGB: Col_1', false);
Map.addLayer(medianAnual_2, rgbMedian,'RGB: Col_2', false);
//Map.addLayer(col_2.map(SAVI).median().clip(municipio),visNDVI,'SAVI');
Map.addLayer(col_2.map(NDFI).median().clip(UMF),visNDFI,'NDFI_depois', false);
Map.addLayer(col_1.map(NDFI).median().clip(UMF),visNDFI,'NDFI_antes',false);
Map.addLayer(col_2.map(DETEX).median().clip(UMF),{'bands': 'DETEX', min:0,max:150, gamma:1},'DETEX_depois');
//Map.addLayer(col_1.map(DETEX).median().clip(UMF),{'bands': 'DETEX', min:1, max:1500, gamma:1},'DETEX_antes',false);
