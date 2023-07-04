osmTiles = function(name, xyz, suffix) {
  result = c(
    osm = "http://tile.openstreetmap.org",
    'osm-admin' = 'http://korona.geog.uni-heidelberg.de/tiles/adminb',
    'osm-roads-grey' = 'http://korona.geog.uni-heidelberg.de/tiles/roadsg/',
    'osm-roads' = 'http://korona.geog.uni-heidelberg.de/tiles/roads',
    'osm-semitransparent' = 'http://korona.geog.uni-heidelberg.de/tiles/hybrid/',
    "osm-no-labels"="http://c.tiles.wmflabs.org/osm-no-labels/",
    "osm-de"="http://c.tile.openstreetmap.de/tiles/osmde/",
    "osm-ru" = "http://a.tiles.wmflabs.org/osm-multilingual/ru,_/",
    "osm-transport"="http://tile2.opencyclemap.org/transport/",
    "stamen-toner" = "https://stamen-tiles-d.a.ssl.fastly.net/toner/",
    "stamen-watercolor" = "https://tiles.stadiamaps.com/styles/stamen_watercolor/",
    "bw-mapnik"="http://b.tiles.wmflabs.org/bw-mapnik2/",
#			mapquest="http://otile1.mqcdn.com/tiles/1.0.0/osm/",
#			"mapquest-sat"="http://otile1.mqcdn.com/tiles/1.0.0/sat",
#      "mapquest-labels"='http://otile3.mqcdn.com/tiles/1.0.0/hyb/',
    'osm-cyclemap' = 'http://a.tile.opencyclemap.org/cycle/',
    'osm-seamap' = 'http://tiles.openseamap.org/seamark/',
    'osm-fr' = 'http://a.tile.openstreetmap.fr/osmfr/',
    'landscape'="http://tile.opencyclemap.org/landscape/",
    'hyda' = 'http://c.tile.openstreetmap.se/hydda/full/',
    'hyda-base' = 'http://c.tile.openstreetmap.se/hydda/base/',
    'hyda-roads' = 'http://c.tile.openstreetmap.se/hydda/roads_and_labels/',
    "opentopomap" = "http://opentopomap.org/",
    "maptoolkit"="http://tile2.maptoolkit.net/terrain/",
    waze="https://worldtiles2.waze.com/tiles/",
    'waze-us'='https://livemap-tiles2.waze.com/tiles/',
    humanitarian="http://a.tile.openstreetmap.fr/hot/",
    cartodb='http://c.basemaps.cartocdn.com/light_all/',
    'cartodb-dark'='http://c.basemaps.cartocdn.com/dark_all/',
#  historical='http://www.openhistoricalmap.org/ohm_tiles/',
    nrcan = 
      'http://geoappext.nrcan.gc.ca/arcgis/rest/services/BaseMaps/CBMT_CBCT_GEOM_3857/MapServer/tile/',
    'nrcan-text' = 
      'http://geoappext.nrcan.gc.ca/arcgis/rest/services/BaseMaps/CBMT_TXT_3857/MapServer/tile/',
    'nrcan-text-fr' = 
      'http://geoappext.nrcan.gc.ca/arcgis/rest/services/BaseMaps/CBCT_TXT_3857/MapServer/tile/',
    spinal = 'http://c.tile.thunderforest.com/spinal-map/',
    neighbourhood = 'https://a.tile.thunderforest.com/neighbourhood/',
    pioneer = 'https://b.tile.thunderforest.com/pioneer/',
    'mobile-atlas'='https://b.tile.thunderforest.com/mobile-atlas/',
    wikimedia = 'https://maps.wikimedia.org/osm-intl/',
    'sputnik' = 'http://tiles.maps.sputnik.ru/'
  )

  # toronto
  #https://gis.toronto.ca/arcgis/rest/services/basemap/cot_topo/MapServer/tile/9/186/142
# https://map.toronto.ca/maps/map.jsp?app=TorontoMaps_v2

#	skobbler="http://tiles3.skobbler.net/osm_tiles2/",	
#		"osm2world"="http://tiles.osm2world.org/osm/pngtiles/n/",
#		bvg="http://mobil.bvg.de/tiles/",
#	landshaded="http://tiles.openpistemap.org/landshaded/",
#		"osm-retina"="http://tile.geofabrik.de/osm_retina/",
#      'osm-rail' = 'http://a.tiles.openrailwaymap.org/standard/',
# rail is 512 insstead of 256 tiles
#			hill="http://www.toolserver.org/~cmarqu/hill/",
#	eu="http://alpha.map1.eu/tiles/",
# 'esri' = 'http://services.arcgisonline.com/ArcGIS/rest/services/World_Street_Map/MapServer/tile/',
# 'esri-grey' = 'http://services.arcgisonline.com/ArcGIS/rest/services/Canvas/World_Light_Gray_Base/MapServer/tile/',
# 'esri-transport'='http://server.arcgisonline.com/ArcGIS/rest/services/Reference/World_Transportation/MapServer/tile/',
# 'esri-topo' = 'http://services.arcgisonline.com/ArcGIS/rest/services/World_Topo_Map/MapServer/tile/'
  
  
# http://server.arcgisonline.com/arcgis/rest/services/Ocean/World_Ocean_Base/MapServer/tile/2/1/1.jpg
  
  # language labels don't appear to be working
  languages = c("en","fr","de", "it","es","ru")
  toadd =	paste("http://a.www.toolserver.org/tiles/osm-labels-", languages,"/", sep="")
  names(toadd) = paste("osm-labels-", languages, sep="")
#	result = c(result, toadd)
  
  
  result = c(result, toadd)
  
  
  
  if(!missing(name)) {
    if(all(name %in% names(result), na.rm=TRUE)) {
      result = result[name]
    } else {
      result = name
    }
  }
  if(!missing(xyz))
    attributes(result)$tileNames = xyz	
  if(!missing(suffix))
    attributes(result)$suffix = suffix	
  
  result
  
}



if(F){
 x = vect(as.matrix(expand.grid(seq(-5e6,-1e6,len=100), seq(-3e6,3e6,len=100))), crs='EPSG:3573')
 xLL = project(x, crsLL)
 zoom=4
 theExt = ext(-6e6,6e6,-6e6,6e6)
 xx = getTilesFromPoints(
  outraster=rast(theExt, res = (xmax(theExt)-xmin(theExt))/1200, crs=crs('EPSG:3031')),
  path="https://stamen-tiles-d.a.ssl.fastly.net/toner/",  
#  path="https://tiles.stadiamaps.com/styles/stamen_watercolor/", suffix='.jpg',
  zoom=3, verbose=TRUE)

}


openmap = function(
  x, 
  zoom, 
  path="http://tile.openstreetmap.org/",
  maxTiles = 9,
  crs=terra::crs(x),
  buffer=0, fact=1,
  verbose=getOption('mapmiscVerbose'),
  cachePath=getOption('mapmiscCachePath')
) {
  


  verbose = max(c(0, verbose))
  
  if(is.null(cachePath)) {
    cachePath = tempdir()
  }
  if(!nchar(cachePath)) {
    cachePath = '.'
  }
  
  alltiles = osmTiles()
  pathOrig = path
  pathisname = gsub("-", ".", pathOrig) %in% 
    gsub("-", ".", names(alltiles))
  path[pathisname] = alltiles[path[pathisname]]
  
  if(length(grep("[[:punct:]]$", path, invert=TRUE)))
    path[ grep("[[:punct:]]$", path, invert=TRUE)] =
      paste(path[ grep("[[:punct:]]$", path, invert=TRUE)], 
        "/", sep="")
  
  if(length(grep("^http[s]*://", path, invert=TRUE)))
    path[ grep("^http[s]*://", path, invert=TRUE)] = 
      paste("http://", 
        path[ grep("^http[s]*://", path, invert=TRUE)], sep="")
  names(path) = pathOrig
  
  
  if(all(class(x) == 'CRS')) {
    # x is a crs object
    # get the ellipse
    crs = x
    toCrop = attributes(x)$ellipse
    x = attributes(x)$regionLL
  } else {
    toCrop = NULL
  }
  


  crsOut=crs

  crsIn = terra::crs(x)

  if(all(is.na(crsIn))) {
    if(is.vector(x)){
      crsIn=crsLL
    } else{
      crsIn = crs	
    }
  }
  

  extMerc = .getExtent(x,crsIn, buffer, crsMerc)
  extMerc = .cropExtent(extMerc, extentMerc)

  if(missing(zoom)) {
    zoom = 1
    while(nTilesMerc(extMerc, zoom) <= maxTiles & zoom <= 18) {
      zoom = zoom + 1
    }
    zoom = min(c(18,max(c(1, zoom-1))))
  }
  if(verbose) cat("zoom is ", zoom, ", ", nTilesMerc(extMerc, zoom), "tiles\n")
  

  # create out raster
  # find average area of pixels in downloaded tiles
  rasterSphere = terra::crop(.getRasterMerc(zoom), extMerc)
  cellSizeMerc = mean(terra::values(terra::cellSize(rasterSphere)))/(256^2)
    # each tile is 256 by 256
  
  # get cell size if testX cells in x direction of output raster
      testX = 100
      testRast = rast(terra::ext(x), res = (terra::xmax(x) - terra::xmin(x))/testX, crs = crsOut)
      testArea = quantile(unlist(terra::spatSample(terra::cellSize(testRast), size=200)), prob=0.5)
      # find number of cells so pixel area matches dowloaded tiles pixel area
      areaRatio = testArea/cellSizeMerc
      newX = fact*round(testX*sqrt(areaRatio))

      outraster = rast(terra::ext(x), res = (terra::xmax(x) - terra::xmin(x))/newX, crs = crsOut)



  
  Dpath = names(path)[1]
  Durl = path[1]

    if(verbose){
      cat(Dpath, '\n')
      cat(Durl, '\n')
    }
    
    if(length(grep(
        'nrcan\\.gc\\.ca|gov\\.bc\\.ca', Durl))
      ){
      suffix = ''
      tileNames = 'zyx'
    } else if(
      length(grep(
          '[.]arcgisonline[.]com', Durl
        ))) {
      suffix='.jpg'
      tileNames = 'zyx'
    } else if(
      length(grep(
          'stamen.watercolor', Durl
        ))) {
      suffix='.jpg'
      tileNames = 'zyx'
    } else if(
      length(grep(
          'heidelberg.de/tiles/(hybrid|adminb|roadsg|roads)/?$', 
          Durl)) |
      length(grep(
          '&$',Durl))
      ){
      tileNames = 'xyz='
      suffix = ''
    } else {
      suffix = '.png'
      tileNames = 'zxy'
    }
    
    if(length(attributes(pathOrig)$tileNames))
      tileNames = attributes(pathOrig)$tileNames
    if(length(attributes(pathOrig)$suffix))
      suffix = attributes(pathOrig)$suffix
    



    result = try(
      getTiles(outraster, 
        zoom=zoom,
        path=Durl,
        verbose=verbose,
        suffix=suffix,
        tileNames = tileNames,
        cachePath = cachePath),
      silent=!verbose)
    
    if(any(class(result)=="try-error")){
      message(paste(Durl, "not accessible"))
    # create an empty raster
    result = outraster
    attributes(result)$openmap = list(
      tiles=NA,
      message=result,
      path=path,
      pathOrig=pathOrig,
      zoom=zoom
    )
  } else {
    attributes(result)$openmap = list(
      path=path,
      pathOrig=pathOrig,
      zoom=zoom
    )
  }
   


  result
}

