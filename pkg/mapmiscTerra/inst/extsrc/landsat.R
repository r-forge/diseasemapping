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
  tilesCrop2 = spTransform(tilesCrop, projection(xExtent))
  tilesCrop3 = raster::crop(tilesCrop2, xExtent)
  
  tiles = tiles[tiles$PR %in% tilesCrop3$PR, ]
  
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
  
  attributes(tilesProjCrop)$landsat = xExtent
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
  
  if(!any(file.size(localFile)>0, na.rm=TRUE)) {
    download.file(url, localFile, mode='wb')
  }
  landsatScenes = read.table(
    gzfile(localFile), sep=',', stringsAsFactors=FALSE,
    header=TRUE
  )
  
  save(landsatScenes, file=rdataFile, compress='xz')
  
  landsatScenes
}

getLandsatCoverage = function(
  path = c(
    options()$mapmiscCachePath, 
    tempdir())[1],
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
  path = c(
    options()$mapmiscCachePath, 
    tempdir())[1],
  mode = 'wb',
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
    if(!any(file.size(scenesLong[D, 'localFile'])>1, na.rm=TRUE) )
      download.file(scenesLong[D,'url'], scenesLong[D, 'localFile'], 
        mode=mode, ...)
  }
  result = tapply(
    X=scenesLong$localFile,
    INDEX=scenesLong[,'entityId', drop=FALSE],
    FUN = function(X)
      raster::stack(sort(X)),
    simplify=FALSE
  )
  lapply(result, function(x) x)
}

landsatNdvi = function(x, 
  suffix = 1:3,
  path = c(
    options()$mapmiscCachePath, 
    tempdir())[1],
  filename = 'ndvi.grd') {
  
  if(basename(filename) == filename) {
    filename = file.path(path, filename)
  }  
  
  suffixString = paste("_(", paste(suffix, collapse='|'), ')$', sep='')
  StilesAll = gsub(suffixString, '', names(x))
  Stiles = unique(StilesAll)
  
  tileMat = list()
  for(D in suffix)
    tileMat[[D]] = outer(
      names(x), 
      paste(Stiles, '_', D, sep=''), 
      '==')
  
  numerMat = tileMat[[suffix[2]]] - tileMat[[suffix[1]]]
  denomMat = tileMat[[suffix[2]]] + tileMat[[suffix[2]]]
  qMat = tileMat[[suffix[3]]]
  
  colnames(numerMat) = colnames(denomMat) = colnames(qMat) =
    Stiles
  
  calc(
    x,
    function(xx) {
      numer = xx %*% numerMat 
      denom = xx %*% denomMat 
      qual = xx %*% qMat 
      qual[qual > 25000 | is.na(numer) | is.na(denom) |
          denom == 0] = NA
      result = numer/denom
      result[is.na(qual)] = NA
      meanResult = mean(result, na.rm=TRUE)
      meanResult[is.nan(meanResult)] = NA
      result = c(result, meanResult, sum(!is.na(result)))
      names(result) = c(colnames(numerMat), 'mean', 'N')
      result
    },
    filename = filename,
    overwrite = file.exists(filename)  
  )
}
#'

