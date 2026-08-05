# Tile layer registry. Stadia (stamen-*) and Thunderforest layers may require
# an API key; append ?apikey=... or /{apikey}/ as required by the provider.
# See openmapAttribution() and the package NEWS.
osmTiles = function(name, xyz, suffix) {
  result = c(
    osm = "https://tile.openstreetmap.org",
    'osm-fr' = 'https://a.tile.openstreetmap.fr/osmfr/',
    "osm-de"="https://tile.openstreetmap.de/",
    "osm-transport"="https://tile.memomaps.de/tilegen/",
    'osm-bw' = "https://tile.openstreetmap.de/tilesbw/osmde/",
    "stamen-toner" = "https://tiles.stadiamaps.com/tiles/stamen_toner/",
    "stamen-watercolor" = "https://watercolormaps.collection.cooperhewitt.org/tile/watercolor/",
    'stamen-terrain' = 'https://tiles.stadiamaps.com/tiles/stamen_terrain/',
    'bvg' = 'https://bvg-gis-c.hafas.de/hafas-tiles/inno2017/2/',
    'esri' = 'https://server.arcgisonline.com/ArcGIS/rest/services/World_Street_Map/MapServer/tile/',
    'esri-satellite' = 'https://server.arcgisonline.com/ArcGIS/rest/services/World_Imagery/MapServer/tile/',
    'esri-natgeo' = 'https://services.arcgisonline.com/ArcGIS/rest/services/NatGeo_World_Map/MapServer/tile/',
    'esri-overlay' = 'https://server.arcgisonline.com/ArcGIS/rest/services/Reference/World_Boundaries_and_Places/MapServer/tile/',
    'esri-topo' = 'https://services.arcgisonline.com/ArcGIS/rest/services/World_Topo_Map/MapServer/tile/',
    kartogiraffe = 'https://tiles.kartogiraffe.de/tiles/map/1',
    'osm-cyclemap' = 'https://c.tile.thunderforest.com/cycle/',
    'osm-seamap' = 'https://tiles.openseamap.org/seamark/',
    # Thunderforest landscape; https has a certificate hostname mismatch
    'landscape'="http://tile.opencyclemap.org/landscape/",
    "opentopomap" = "https://a.tile.opentopomap.org/",
    waze="https://worldtiles3.waze.com/tiles/",
    'waze-us'='https://livemap-tiles1.waze.com/tiles/',
    humanitarian="https://a.tile.openstreetmap.fr/hot/",
    cartodb = 'https://basemaps.cartocdn.com/light_all/',
    'cartodb-nolabels' = 'https://basemaps.cartocdn.com/light_nolabels/',
    'cartodb-dark'='https://basemaps.cartocdn.com/dark_all/',
    'f4' = 'https://tile3.f4map.com/tiles/f4_2d/',
    nrcan = 'https://maps-cartes.services.geo.ca/server2_serveur2/rest/services/BaseMaps/CBMT_CBCT_GEOM_3857/MapServer/tile/',
    'nrcan-text' =
    'https://maps-cartes.services.geo.ca/server2_serveur2/rest/services/BaseMaps/CBMT_TXT_3857/MapServer/tile/',
    'nrcan-text-fr' =
    'https://maps-cartes.services.geo.ca/server2_serveur2/rest/services/BaseMaps/CBCT_TXT_3857/MapServer/tile/',
    spinal = 'https://c.tile.thunderforest.com/spinal-map/',
    neighbourhood = 'https://a.tile.thunderforest.com/neighbourhood/',
    pioneer = 'https://b.tile.thunderforest.com/pioneer/',
    'mobile-atlas'='https://b.tile.thunderforest.com/mobile-atlas/',
    ump = 'https://2.tiles.ump.waw.pl/ump_tiles/'
  )

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



openmap = function(
  x, 
  zoom, 
  path="https://tile.openstreetmap.org/",
  maxTiles = 9,
  crs=ifelse(is.numeric(x), mapmisc::crsLL, terra::crs(x)),
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
  


  if(identical(crs, "")) {
    crs=crsIn=crsOut = crsLL
  } else {
    crsOut=crs
    crsIn = try(terra::crs(x), silent=TRUE)    
    if(identical(crsIn, "") | any(class(crsIn) == 'try-error')) crsIn = crs
  }
  if(is.numeric(x)) {
    x = vect(matrix(x, ncol=2), crs=crsIn)
  }

  if(all(class(x) == 'SpatExtent')) {
    x = rast(extent = x, crs = crsIn)
  }

# get extent of output
# this ignores negatives
  if(terra::is.lonlat(crsIn)) {
    outExtent = terra::extend(terra::ext(x), buffer)
  } else {
    extVect = as.vector(terra::ext(x))
    extVect= extVect + c(-1,1,-1,1) * rep_len(buffer, 4)
    outExtent = terra::ext(extVect)
  }

  if(!identical(as.character(crsIn), as.character(crsOut))) {

    extentTestRast = outExtent
    if(any(diff(as.vector(extentTestRast))[-2] <= 1e-4)) {
      extentTestRast = terra::extend(extentTestRast, 1e-4)
    }
    testRast = rast(extentTestRast, res = (terra::xmax(extentTestRast) - terra::xmin(extentTestRast))/NtestCols, crs = crsIn)
    testPoints = vect(terra::xyFromCell(testRast, 1:terra::ncell(testRast)), crs=terra::crs(testRast))

    testPointsOut = suppressWarnings(terra::project(testPoints, crsOut))
    outExtent= terra::ext(testPointsOut)

  }




  if(terra::is.lonlat(crsOut)) {
    outExtent = terra::intersect(outExtent, terra::unwrap(bboxLL))
  }

  testRast = rast(outExtent, 
    res = (terra::xmax(outExtent) - terra::xmin(outExtent))/NtestCols, crs = crsOut)
  testPoints = suppressWarnings(terra::as.points(testRast))
  testPointsMerc = suppressWarnings(terra::project(testPoints, crsMerc))

  if(missing(zoom)) {
    # get zoom

    zoom = 0
    Ntiles = 0
    while(Ntiles <= maxTiles & zoom <= 18) {
      zoom = zoom + 1
      Ntilesm1 = Ntiles
      Ntiles = length(unique(terra::cellFromXY(
        .getRasterMerc(zoom), terra::crds(testPointsMerc))))
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
  'nrcan\\.gc\\.ca|gov\\.bc\\.ca|maps-cartes\\.services\\.geo\\.ca', Durl))
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
    'stamen.watercolor|watercolormaps[.]collection[.]cooperhewitt[.]org', Durl
  ))) {
  suffix='.jpg'
  tileNames = 'zxy'
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


result = try(
  getTiles(outraster, 
    zoom=zoom,
    path=Durl,
    verbose=verbose,
    suffix=suffix,
    tileNames = tileNames,
    cachePath = cachePath),
  silent=!verbose)

if(inherits(result, 'try-error')){
    # create an empty raster
  result = outraster
  terra::values(outraster) = NA
}

if(all(is.na(terra::values(result)))){
  message(paste(Durl, "not accessible"))
    # create an empty raster
  attributes(result)$openmap = list(
    tiles=NA,
    message=result,
    path=path,
    pathOrig=pathOrig,
    zoom=zoom
  )
  return(result)
}


  # fill in poles
thePoles = as.matrix(expand.grid(seq(-170, 180, len=100), 
  as.vector(outer(c(-1,1),seq(83, 90, len=201)))))
thePoles = vect(thePoles, crs=crsLL, 
  atts = data.frame(pole=c('south','north')[1+(thePoles[,2]>0)]))
thePolesTrans = suppressWarnings(terra::project(thePoles, crs(result)))
thePolesTrans$cell = terra::cellFromXY(result, terra::crds(thePolesTrans))
thePolesTrans = thePolesTrans[!is.na(thePolesTrans$cell)]
thePolesTrans = thePolesTrans[!duplicated(thePolesTrans$cell), ,drop=FALSE]

if(length(thePolesTrans)) {
  thePolesSplit = terra::split(thePolesTrans, thePolesTrans$pole)
  names(thePolesSplit) = unlist(lapply(thePolesSplit, function(xx) xx$pole[1]))


  theBuffer = suppressWarnings(lapply(thePolesSplit, function(xx) {
    unique(as.vector(terra::adjacent(result, xx$cell, directions=16, include=TRUE)))
  } ))



  # replace NA's near north pole with values nearby
  for(Ddirection in names(theBuffer)) {
    if(verbose) cat('poles: ', Ddirection, '\n')
      extractHere = terra::extract(result, y=theBuffer[[Ddirection]])
    cellsToFill = theBuffer[[Ddirection]][is.na(extractHere[,1])]
    toFill = na.omit(extractHere[,1])

    if(length(toFill) & length(cellsToFill)) {
      toFill = as.numeric(names(sort(table(toFill), decreasing=TRUE))[1]) 
      result[cellsToFill] = toFill
    }
  }
}


  # if there's only one layer, and no colortable
  # convert to greyscale
if(all(terra::nlyr(result)== 1) & !any(terra::has.colors(result))){
  theRange = unlist(terra::global(result, range))
  if( all(ceiling(theRange) == floor(theRange), na.rm=TRUE) & all(theRange >= 0) & all(theRange < 256)) {
      # 264 grey scale
    terra::coltab(result) = data.frame(value = seq(0,255), col = grDevices::grey(seq(0,1,len=256)) )
  } 
}

  # if there's transparency convert to color table
if(all(c('red','green','blue','alpha') %in% names(result) & 
  terra::nlyr(result) == 4 & 
  !any(terra::has.colors(result)) )) {
  theColTab = terra::unique(result)
  theColTab = cbind(ID = 1:nrow(theColTab), theColTab)

  valuesOrig = terra::values(result, dataframe=TRUE)
  valuesOrig$idOrig = 1:nrow(valuesOrig)

  newValues = merge(valuesOrig, theColTab)
  newValues = newValues[order(newValues$idOrig), 'ID']

  result2 = rast(result, nlyrs=1)
  terra::values(result2) = newValues
  terra::coltab(result2) = theColTab
  result = result2
}

result = writeRasterMapTiles(result, 
  filename = tempfile(tmpdir=cachePath, fileext='.tif'))

attributes(result)$openmap = list(
  path=path,
  pathOrig=pathOrig,
  zoom=zoom
)

result
}

