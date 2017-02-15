#+ packages
if(FALSE) {
  install.packages("Pmisc", 
    repos='http://r-forge.r-project.org')
  # also needs R.utils
}
#'

#+ getLandsatTilesFunction


getLandsatTiles = function(x,
  tiles = getLandsatCoverage(),
  buffer = 1000,
  crsOut = raster::projection(x)) {
  
  xExtent = 
    raster(extend(extent(x), buffer), 
      ncol=10, nrow=10, crs=projection(x))
  
  xExtentProj = projectExtent(xExtent, projection(tiles))
  
  tilesCrop = raster::crop(
    tiles, xExtentProj
  )
  tiles = tiles[tiles$PR %in% tilesCrop$PR, ]
  
  if(identical(crsOut, 'landsat'))
    crsOut = paste(
      '+proj=lsat +path=',
      tiles$PATH[1],
      ' +lsat=5 +ellps=WGS84 +datum=NAD83', sep='')
  
  if(identical(crsOut, 'omerc'))
    crsOut = mapmisc::omerc(
      tiles, 
      angle=seq(-25,25,by=0.1)
    )
  
  
  tilesProjCrop = spTransform(tiles, crsOut)
  
  tilesProjCrop$pPath = paste('p',  sprintf("%03d", tilesProjCrop$PATH), sep='')
  tilesProjCrop$rRow = paste('r',  sprintf("%03d", tilesProjCrop$ROW), sep='')
  
  tilesProjCrop$pathString = sprintf("%03d", tilesProjCrop$PATH)
  tilesProjCrop$rowString = sprintf("%03d", tilesProjCrop$ROW)
  
  tilesProjCrop$tileString = paste(
    tilesProjCrop$pathString,
    tilesProjCrop$rowString, sep='')
  
  tilesProjCrop$tTile = paste(
    tilesProjCrop$pPath,
    tilesProjCrop$rRow, sep=' ')
  
  tilesProjCrop
}

getLandsatScenes = function(
  timeRange,
  tiles,
  sceneList = getLandsatSceneList()
) {
  
  if(any(names(tiles) == 'tileString')) {
    tiles = tiles$tileString
  }
  
  sceneList$tile = sprintf("%03d%03d", sceneList$path, sceneList$row)
  toKeep = intersect(
    which(sceneList$tile %in% tiles),
    which(sceneList$acquisitionDate >= min(timeRange) &
        sceneList$acquisitionDate <= max(timeRange))
  )
  sceneList[toKeep,]  
}

getLandsatSceneList = function(
  path = gsub("[[:alnum:]]+$", "landsat", tempdir()),
  url = "http://landsat-pds.s3.amazonaws.com/scene_list.gz") {
  
  localFile = file.path(path, basename(url))
  outFile = gsub("[.]gz$", "", localFile, ignore.case = TRUE)
  rdataFile = file.path(path, 'landsatScenes.RData')
  
  if(file.exists(rdataFile)) {
    rdataEnv = new.env()
    load(rdataFile, envir=rdataEnv)
    if('landsatScenes' %in% ls(envir=rdataEnv)) {
      return(get('landsatScenes', envir=rdataEnv))
    }
  }
  
  if(!file.exists(outFile)) {
    if(!file.exists(localFile)) {
      download.file(url, localFile, mode='wb')
    }
    R.utils::gunzip(localFile, 
      remove=FALSE,
      overwrite = file.exists(outFile))
  }
  landsatScenes = read.table(
    outFile, sep=',', stringsAsFactors=FALSE,
    header=TRUE
  )
  
  save(landsatScenes, file=rdataFile, compress='xz')
  landsatScenes
}

getLandsatCoverage = function(
  path = gsub("[[:alnum:]]+$", "landsat", tempdir()),
  url = 'https://landsat.usgs.gov/sites/default/files/documents/wrs2_descending.zip'
) {
  
  localFile = file.path(path, basename(url))
  rdataFile = file.path(path, 'landsatCover.RData')
  if(file.exists(rdataFile)) {
    rdataEnv = new.env()
    load(rdataFile, envir=rdataEnv)
    if('landsatCover' %in% ls(envir=rdataEnv)) {
      return(get('landsatCover', envir=rdataEnv))
    }
  }
  if(!file.exists(localFile)) {
    download.file(url, localFile, mode='wb')
  }
  
  sfile = file.path(path,
    grep("shp$", unzip(localFile, list=TRUE)$Name, value=TRUE)
  )
  
  if(!file.exists(sfile))
    utils::unzip(localFile, exdir=path)
  landsatCoverLL = raster::shapefile(sfile)
  if(requireNamespace('rgdal')) {
    landsatCover = spTransform(
      raster::crop(landsatCoverLL, extent(-180,180,-88,88)),
#        "+proj=eck2 +ellps=WGS84 +datum=NAD83"
      "+proj=sinu +ellps=WGS84 +datum=NAD83"
    )
  } else {
    landsatCover = landsatCoverLL
  }
  save(landsatCover, file=rdataFile, compress='xz')
  landsatCover
}
downloadLandsatTiles = function(
  scenes, band = 1:11,
  path = gsub("[[:alnum:]]+$", "landsat", tempdir()),
  ...
) {
  
  scenesLong = scenes[
    rep(1:nrow(scenes), length(band)),
  ]
  scenesLong$band = rep(band, rep(nrow(scenes), length(band)))
  scenesLong$url = paste(gsub("index[.]html$", "", 
      scenesLong$download_url), 
    scenesLong$entityId, '_B', scenesLong$band, '.TIF', sep='')
  scenesLong$localFile = file.path(path, basename(scenesLong$url))
  
  if(!dir.exists(path)) dir.create(path)
  for(D in 1:nrow(scenesLong)) {
    if(!file.exists(scenesLong[D, 'localFile']))
      download.file(scenesLong[D,'url'], scenesLong[D, 'localFile'], 
        mode='wb', ...)
  }
  result = tapply(
    scenesLong$localFile,
    scenesLong$entityId,
    function(x, ...)
      raster::stack(sort(x), ...)
  )
  result
}

landsatNdvi = function(x, 
  ndviFile = gsub(
    '_B[[:digit:]][.]TIF$', '_NDVI.grd',
    filename(x[[1]]), ignore.case=TRUE)
) {
  
  fourLayer = grep("_B4$", names(x), value=TRUE)
  fiveLayer = grep("_B5$", names(x), value=TRUE)
  # https://landsat.usgs.gov/qualityband
  qaLayer = grep("_BQA$", names(x), value=TRUE)
  
  result = overlay(
    x[[fourLayer]], x[[fiveLayer]], x[[qaLayer]], 
    fun = function(a,b,c) {
      isNA = c(1, NA)[1+(c > 50000 | a==1 | b==1)]
      as.numeric(isNA*(b-a)/(b+a))
    },
    filename=ndviFile, overwrite=file.exists(ndviFile))
  names(result) = gsub("_B[[:digit:]]$", "_NDVI", fourLayer)
  
  result
  
}

#'

#+ setup
dataDir = '/store/patrick/landsat'
if(!dir.exists(dataDir)) dir.create(dataDir)
#'

#+ downloadLandsatCover
library(raster)
landsatCover = getLandsatCoverage(dataDir)
landsatScenes = getLandsatSceneList(dataDir)
#'


#+ landsatToronto
toLL = mapmisc::geocode("Toronto, ON")

landsatTiles = getLandsatTiles(
  x=toLL, 
  tiles = landsatCover, 
  buffer=0.25, crsOut='omerc')

torontoTiles = mapmisc::tonerToTrans(
  mapmisc::openmap(landsatTiles, path='stamen-toner')
)
#'

#+ torontoLandastPlot
Stiles = seq(1, len=length(landsatTiles))
Scol = mapmisc::col2html(
  RColorBrewer::brewer.pal(
    length(landsatTiles), 
    'Dark2'
  ),
  0.85
)
mapmisc::map.new(landsatTiles, buffer=-2*1000)
plot(torontoTiles, add=TRUE)
plot(landsatTiles, add=TRUE, 
  border=Scol,
  col = Scol, 
  density = 12,
  angle = 180*Stiles/length(Stiles),
  lwd=3
)
mapmisc::legendBreaks("bottomleft", 
  breaks= landsatTiles$tTile, 
  col = Scol)    
#'

#+ selectTiles

scenesFeb = getLandsatScenes(
  timeRange = ISOdate(
    2014, 2, c(1,28), tz='utc'
  ),
  tiles = landsatTiles,
  sceneList = landsatScenes
)
#'

#+ downloadTiles
landsatFeb = downloadLandsatTiles(
  scenes=scenesFeb, band = c(4,5,'QA'),
  path = dataDir
)  
#'

#+ plotBand4
torontoTiles2 = mapmisc::tonerToTrans(
  mapmisc::openmap(
    landsatFeb[[1]], 
    path='stamen-toner', 
    fact = 1.5)
)

febCol = mapmisc::colourScale(
  landsatFeb[[2]][[1]], breaks=12, dec=-3, col='Spectral',
  style= 'quantile'
)

plot(landsatFeb[[2]][[1]], col=febCol$col, breaks=febCol$breaks, legend=FALSE)
plot(torontoTiles2, add=TRUE)
mapmisc::legendBreaks("right", febCol)


mapmisc::map.new(good, buffer=-c(20,80,170,10)*1000)
plot(good, col=febCol$col, breaks=febCol$breaks, legend=FALSE, add=TRUE,
  maxpixels=2000000)
plot(torontoTiles2, add=TRUE)
mapmisc::legendBreaks("right", febCol)

#'

#+ calcNdvi
ndviFeb = lapply(landsatFeb, landsatNdvi)
#'

#+ plotNdviFeb
myCol = mapmisc::colourScale(
  ndviFeb[[2]], 
  breaks=c(-1, -0.1, 0, 0.01, 0.02, 0.05, 0.1, 0.2,0.4, 1), 
  col='Spectral',
  style= 'fixed'
  )

mapmisc::map.new(ndviFeb[[2]], buffer=-c(20,80,170,10)*1000)
plot(ndviFeb[[2]], col=myCol$col, breaks=myCol$breaks, legend=FALSE, add=TRUE,
  maxpixels=3000000)
plot(torontoTiles2, add=TRUE)
mapmisc::legendBreaks("right", myCol)
  

#'