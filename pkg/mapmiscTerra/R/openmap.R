# stamen https://maps.stamen.com/stadia-partnership/

# Name  Bucket  Prefix  Extension
#Terrain tile.stamen.com /terrain  png
#Terrain Background  tile.stamen.com /terrain-background png
#Terrain Labels  tile.stamen.com /terrain-labels png
#Terrain Lines tile.stamen.com /terrain-lines  png
#Toner long-term.cache.maps.stamen.com /toner  png
#Toner Background  long-term.cache.maps.stamen.com /toner-background png
#Toner Labels  long-term.cache.maps.stamen.com /toner-labels png
#Toner Lines long-term.cache.maps.stamen.com /toner-lines  png
#Toner Lite  long-term.cache.maps.stamen.com /toner-lite png
#Watercolor  long-term.cache.maps.stamen.com /watercolor jpg
#Each tile is accessible at an S3 URL with the following format:

#s3://{bucket}{prefix}/{z}/{x}/{y}.{extension}
#For example, to get the watercolor tile for zoom 2, x 3, and y 1 you would use:

#s3://long-term.cache.maps.stamen.com/watercolor/2/3/1.jpg


# https://mc.bbbike.org/mc/?num=2&mt0=mapnik&mt1=maptiler_streets
osmTiles = function(name, xyz, suffix) {
  result = c(
    osm = "http://tile.openstreetmap.org",
    'osm-admin' = 'http://korona.geog.uni-heidelberg.de/tiles/adminb',
    'osm-roads-grey' = 'http://korona.geog.uni-heidelberg.de/tiles/roadsg/',
    'osm-roads' = 'http://korona.geog.uni-heidelberg.de/tiles/roads',
    'osm-semitransparent' = 'http://korona.geog.uni-heidelberg.de/tiles/hybrid/',
#    "osm-no-labels"="http://c.tiles.wmflabs.org/osm-no-labels/",
    "osm-de"="http://c.tile.openstreetmap.de/tiles/osmde/",
    "osm-ru" = "http://a.tiles.wmflabs.org/osm-multilingual/ru,_/",
    "osm-transport"="http://tile2.opencyclemap.org/transport/",
    "stamen-toner" = "https://tiles.stadiamaps.com/tiles/stamen_toner/",
    "stamen-watercolor" = "https://watercolormaps.collection.cooperhewitt.org/tile/watercolor/",
    'stamen-terrain' = 'https://stamen-tiles-c.a.ssl.fastly.net/terrain/',
    "bw-mapnik"="http://b.tiles.wmflabs.org/bw-mapnik2/",
    'bvg' = 'https://bvg-gis-c.hafas.de/hafas-tiles/inno2017/',
    'esri' = 'https://server.arcgisonline.com/ArcGIS/rest/services/World_Street_Map/MapServer/tile/',
    'esri-satellite' = 'https://server.arcgisonline.com/ArcGIS/rest/services/World_Imagery/MapServer/tile/',
    'esri-natgeo' = 'https://services.arcgisonline.com/ArcGIS/rest/services/NatGeo_World_Map/MapServer/tile/',
    'esri-overlay' = 'https://server.arcgisonline.com/ArcGIS/rest/services/Reference/World_Reference_Overlay/MapServer/tile/',
    'esri-topo' = 'https://services.arcgisonline.com/ArcGIS/rest/services/World_Topo_Map/MapServer/tile/',
    'komoot' = 'https://a.tile.hosted.thunderforest.com/komoot-2/',
#    soviet = 'https://y.tile.bbbike.org/cgi-bin/tapp/tilecache.py/1.0.0/topomapper_v2/',
    'osm-cyclemap' = 'http://a.tile.opencyclemap.org/cycle/',
    'osm-seamap' = 'http://tiles.openseamap.org/seamark/',
    'osm-fr' = 'http://a.tile.openstreetmap.fr/osmfr/',
    'landscape'="http://tile.opencyclemap.org/landscape/",
    'hyda' = 'http://c.tile.openstreetmap.se/hydda/full/',
    'hyda-base' = 'http://c.tile.openstreetmap.se/hydda/base/',
    'hyda-roads' = 'http://c.tile.openstreetmap.se/hydda/roads_and_labels/',
    "opentopomap" = "https://a.tile.opentopomap.org/",
    "maptoolkit"="https://rtc-cdn.maptoolkit.net/rtc/toursprung-terrain/",
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



if(FALSE){
 x = vect(as.matrix(expand.grid(seq(-5e6,-1e6,len=100), seq(-3e6,3e6,len=100))), crs='EPSG:3573')
 xLL = project(x, crsLL)
 zoom=4
 theExt = ext(-6e6,6e6,-6e6,6e6)
 xx = openmap(
  x=rast(theExt, res = (xmax(theExt)-xmin(theExt))/1200, crs=crs('EPSG:3031')),
  path="https://stamen-tiles-d.a.ssl.fastly.net/toner/",  
#  path="https://tiles.stadiamaps.com/styles/stamen_watercolor/", suffix='.jpg',
  zoom=2, verbose=TRUE, fact=2)

}


openmap = function(
  x, 
  zoom, 
  path="http://tile.openstreetmap.org/",
  maxTiles = 9,
  crs=terra::crs(x),
  buffer=0, fact=1,
  verbose=getOption('mapmiscVerbose'),
  cachePath=getOption('mapmiscCachePath'), 
  suffix=NULL
) {



  verbose = max(c(0, verbose))
  

  NtestCols = 100

  
  if(!is.null(attributes(x)$ellipse) ) {
    # to do: check for ellipses
    # x is a crs object
    # get the ellipse
    crs = x
    toCrop = attributes(x)$ellipse
    x = attributes(x)$ellipse
  } else {
    toCrop = NULL
  }
  

  if(is.numeric(x)) x = vect(matrix(x, ncol=2), crs=crs)

  if(all(class(x) == 'SpatExtent')) x = rast(extent = x, crs = crs)

  if(identical(crs, "")) {
      crs=crsIn=crsOut = crsLL
  } else {
    crsOut=crs
    crsIn = terra::crs(x)    
  }
  

# get extent of output

# get output raster
  extentTestRast = terra::extend(terra::ext(x), buffer)
  if(any(diff(as.vector(extentTestRast))[-2] <= 1e-4)) {
      extentTestRast = terra::extend(extentTestRast, 1e-4)
  }
  testRast = rast(extentTestRast, res = (terra::xmax(extentTestRast) - terra::xmin(extentTestRast))/NtestCols, crs = crsIn)
  testPoints = vect(terra::xyFromCell(testRast, 1:terra::ncell(testRast)), crs=terra::crs(testRast))

  testPointsMerc = project(testPoints, crsMerc)
  testPointsOut = terra::project(testPoints, crsOut)
  outExtent= terra::ext(testPointsOut)

# buffer
  if(terra::is.lonlat(crsOut) & any(buffer > 90)) {
          # assume buffer is in km
          # transform to merc , buffer, transform back
          outExtentMerc = terra::extend(terra::ext(testPointsMerc), buffer)
          outPointsMerc = vect(matrix(as.vector(outExtentMerc), ncol=2), crsMerc)
          outPointsLL = project(outPointsMerc, crsOut)
          outExtent = terra::ext(outPointsLL)

    } else {
          outExtent = terra::extend(outExtent, buffer)
    }
    if(terra::is.lonlat(crsOut)) {
      outExtent = terra::intersect(outExtent, terra::unwrap(bboxLL))
    }

  testRast = rast(outExtent, res = (terra::xmax(outExtent) - terra::xmin(outExtent))/NtestCols, crs = crsOut)
  testPoints = vect(terra::xyFromCell(testRast, 1:terra::ncell(testRast)), crs=terra::crs(testRast))
  testPointsMerc = project(testPoints, crsMerc)

  if(missing(zoom)) {


    # get zoom

    zoom = 0
    Ntiles = 0
    while(Ntiles <= maxTiles & zoom <= 18) {
      zoom = zoom + 1
      Ntilesm1 = Ntiles
      Ntiles = length(unique(terra::cellFromXY(.getRasterMerc(zoom), terra::crds(testPointsMerc))))
    }
    zoom = zoom - 1
    if(verbose) cat("zoom is ", zoom, ", ", Ntilesm1, "tiles\n")
  }



  # create out raster
  # find average area of pixels in downloaded tiles

    mercHere = .getRasterMerc(zoom)
   if(identical(crsOut, crsMerc)) {
        # output crs is mercator, return tiles as-is
        outraster = terra::crop(mercHere, testRast, snap='out')
        outraster = terra::disagg(outraster, 256)

    } else {
    theTable = as.data.frame(table(terra::cellFromXY(mercHere, terra::crds(testPointsMerc))))
    theTable$cell = as.numeric(as.character(theTable[,1]))
 
    mercHere = terra::crop(mercHere, 
      terra::ext(
        rep(terra::xyFromCell(mercHere, theTable[which.max(theTable$Freq), 'cell']), each=2) + 
        0.6*rep(terra::res(mercHere), each=2)*c(-1,1,-1,1)
      )
    )
    # each tile is 256 x 256
    mercHere = terra::disagg(mercHere, 256)

    # width of the cell with the most test points in it
    cellWidthMerc = quantile(terra::values(terra::cellSize(mercHere, unit='m')), 0.5, na.rm=TRUE)


    areaRast = suppressWarnings(terra::cellSize(testRast, unit='m'))
    if(terra::ncell(areaRast) < 1e4) {
      toQuantile = terra::values(areaRast)
    } else {
      toQuantile = unlist(terra::spatSample(areaRast,  size=min(c(terra::ncell(areaRast), 1e4))))
    }
    toQuantile = toQuantile[toQuantile > 0]
    cellWidthRast = quantile(toQuantile, prob=0.5, na.rm=TRUE)


    areaRatio = cellWidthRast/cellWidthMerc

    newNumberOfCells = fact*NtestCols*sqrt(areaRatio)

    outraster = rast(outExtent, res = (terra::xmax(outExtent) - terra::xmin(outExtent))/newNumberOfCells, crs = crsOut)
  } # end not merc

# cache


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
  
  
  Dpath = names(path)[1]
  Durl = path[1]

  suffixOrig = suffix

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
      '[.]arcgisonline[.]com|bbbike', Durl
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
  if(!is.null(suffixOrig))
    suffix = suffixOrig


#  stuff <<- list(outraster, zoom, Durl, verbose, suffix, tileNames, cachePath)
  # outraster = stuff[[1]];zoom=stuff[[2]]; path = stuff[[3]]; verbose=stuff[[4]]; suffix=stuff[[5]]; tileNames = stuff[[6]]; cachePath = stuff[[7]]
  result = try(
    getTiles(outraster, 
      zoom=zoom,
      path=Durl,
      verbose=verbose,
      suffix=suffix,
      tileNames = tileNames,
      cachePath = cachePath),
    silent=!verbose)

  # fill in poles
  # TO DO: doesn't work with omerc
  thePoles = as.matrix(expand.grid(seq(-170, 180, by=10), as.vector(outer(c(-1,1),c(84.5, 85.5)))))
  thePolesTrans = project(vect(thePoles, crs=crsLL, atts = data.frame(pole=c('south','north')[1+(thePoles[,2]>0)])), crs(result))
  thePolesSplit = terra::split(thePolesTrans, thePolesTrans$pole)
  names(thePolesSplit) = unlist(lapply(thePolesSplit, function(xx) xx$pole[1]))
  theHull = suppressWarnings(lapply(thePolesSplit, terra::convHull))
  theHull = theHull[unlist(lapply(theHull, length))>1]

  # replace NA's near north pole with values nearby
  if(!is.null(theHull$north)) {
    toFill = terra::extract(result, thePolesTrans[thePolesTrans$pole=='north'], ID=FALSE)
    toFill = toFill[min(which(!is.na(toFill[,1]))),,drop=FALSE]

    result = terra::mask(result, theHull$north, updatevalue = (toFill), inverse=TRUE)
  }
  if(!is.null(theHull$south)) {
    toFill = terra::extract(result, thePolesTrans[thePolesTrans$pole=='south'], ID=FALSE)
    toFill = toFill[min(which(!is.na(toFill[,1]))),,drop=FALSE]

    result= terra::mask(result, theHull$south, updatevalue = unlist(toFill), inverse=TRUE)
  }

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

