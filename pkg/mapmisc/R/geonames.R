GNcities = function(north, east, south, west, lang = "en", 
  maxRows = 10, buffer=0) {
  
  if(length(options()$mapmiscVerbose)) { 
    verbose = options()$mapmiscVerbose
  } else {
    verbose=FALSE
  }
  
  cachePath=options()$mapmiscCachePath
  
  if(is.null(cachePath)) {
    cachePath = tempdir()
  }
  if(!nchar(cachePath)) {
    cachePath = '.'
  }
  
  cachePath = file.path(cachePath,'GNcities')
  dir.create(cachePath,recursive=TRUE,showWarnings=FALSE)
  
  fourCoords=FALSE
  if(is.numeric(north))
    if(length(north)== 1)
      fourCoords = TRUE
  
  theproj = projection(north)
  if(!fourCoords) {
    extLL = .getExtent(north, extend=buffer, crsOut = crsLL)
    
    east = xmax(extLL)
    west = xmin(extLL)
    south = ymin(extLL)
    north= ymax(extLL)
    
  }
  
  cacheFile = file.path(
    cachePath,  
    paste(
      make.names(toString(c(north, east, south,west,lang,maxRows))), 
      '.Rdata', sep='')
  )
  
  if(file.exists(cacheFile)) {
    print(1)
    if(verbose) {
      message("found in cache", cachePath)
    }
    load(cacheFile)
  }	else if(requireNamespace("geonames", quietly = TRUE)) { 
    if(verbose)	message("not found in cache, retrieving")
    result = geonames::GNcities(north=north,east=east,
      south=south,west=west,lang,maxRows)
    if(verbose)
      message("caching in", cachePath)
    save(result, file=cacheFile)
  } else {
    warning("install the geonames package to use GNcities")
    result = NULL
  }
  result = SpatialPointsDataFrame(cbind(
      as.numeric(result[,'lng']),
      as.numeric(result[,'lat'])
    ), data=result, 
    proj4string=mapmisc::crsLL)
  
  if(is.na(theproj)) {
    if(requireNamespace('rgdal', quietly=TRUE ))
      result = spTransform(result, CRSobj=CRS(theproj))
  }
  
  result
}

GNsearch = function(..., crs=crsLL) {
  
  
  if(requireNamespace("geonames", quietly = TRUE)) {
    
    theDots = list(...)
    isVector = unlist(lapply(theDots, length))
    isVector = isVector[isVector > 1]
    
    
    
    if(length(isVector)) {
      result = mapply(
        geonames::GNsearch,
        ...,
        SIMPLIFY=FALSE
      )
      result = do.call(rbind, result)    
      result = as.data.frame(result)
    } else {
      result=geonames::GNsearch(...)
    }
    
    if(all(c("lat","lng") %in% names(result))){
      coords = as.matrix(result[,c("lng","lat"),drop=FALSE])
      mode(coords) = 'numeric'
      
      result$population = as.numeric(result$population)
      
      result = SpatialPointsDataFrame(
        coords,
        data=result, 
        proj4string=crsLL)
      if(requireNamespace('rgdal', quietly=TRUE ))
        result = spTransform(result, CRSobj=crs(crs))
    }
    
  } else {
    warning("install the geonames package to use GNsearch")
    result = NULL
    
  }
  result
}


geocode = function(x, oneRecord=FALSE, extent=NULL, progress='', ...) {
  
  if(length(options()$mapmiscVerbose)) { 
    verbose = options()$mapmiscVerbose
  } else {
    verbose=FALSE
  }
  cachePath=options()$mapmiscCachePath
  if(is.null(cachePath)) {
    cachePath = tempdir()
  }
  if(!nchar(cachePath)) {
    cachePath = '.'
  }
  cachePath = file.path(cachePath,'geocode')
  dir.create(cachePath,recursive=TRUE,showWarnings=FALSE)
  
  theproj = projection(extent)
  extLL = .getExtent(extent, crsOut = crsLL)
  
  cacheFile = file.path(
    cachePath,  
    paste(
      make.names(toString(c(x,as.vector(extLL), oneRecord))), 
      '.Rdata', sep='')
  )
  
  if(file.exists(cacheFile)) {
    if(verbose)
      message("found in cache ", cachePath)
    load(cacheFile)
  }	else if(requireNamespace("dismo", quietly = TRUE)) {
    if(verbose)
      message("not found in cache, retrieving")
    result = dismo::geocode(
      x, oneRecord=oneRecord, 
      extent=extLL, progress=progress, ...)		
    if(verbose)
      message("caching in ", cachePath)
    
    saveRes = try(save(result, file=cacheFile), silent=TRUE)
    if(verbose & class(saveRes) == 'try-error') {
      message("cached file couldnt be saved")
    }
  } else {
    warning("install the dismo package to use geocode")
    result = NULL
  }
  
  if(is.data.frame(result)) {
    result$name = gsub(", [[:print:]]+$", "", 
      as.character(result$interpretedPlace))
    resultCoords = as.matrix(result[,c('longitude','latitude')])
    result = SpatialPointsDataFrame(
      resultCoords,
      data = result,
      proj4string = mapmisc::crsLL
    )
  }
  result
  
}
