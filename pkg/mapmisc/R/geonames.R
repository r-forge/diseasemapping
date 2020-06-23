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
  dirCreateMapmisc(cachePath,recursive=TRUE,showWarnings=FALSE)
  
  fourCoords=FALSE
  if(is.numeric(north))
    if(length(north)== 1)
      fourCoords = TRUE
  
  theproj = crs(north)
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
    if(verbose) {
      message("found in cache", cachePath)
    }
    load(cacheFile)
  }	else if(requireNamespace("geonames", quietly = TRUE)) { 
     if(verbose)	message("not found in cache, retrieving")
    result = mapmiscGNcities(north=north,east=east,
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


geocode = function(x, 
  extent,
  lang = gsub("(_|[:]).*", "", Sys.getenv('LANGUAGE'))
  ) {
#  x = paste(c('ottawa','nain','winnipeg'), 'canada', sep=',')

  if(length(getOption('mapmiscVerbose'))) { 
    verbose = getOption('mapmiscVerbose')
  } else {
    verbose=FALSE
  }

  cachePath=getOption('mapmiscCachePath')
  if(is.null(cachePath)) {
    cachePath = tempdir()
  }
  if(!nchar(cachePath)) {
    cachePath = '.'
  }
  cachePath = file.path(cachePath,'geocode')
  dirCreateMapmisc(cachePath,recursive=TRUE,showWarnings=FALSE)

  if(!nchar(lang)) {
    lang = 'en'
  }
  langString = paste0(paste(c('name', lang), collapse=':'), ": ")

  xDf = data.frame(
  orig = x, 
  url=paste0(
    'https://nominatim.openstreetmap.org/search/',
    gsub("[[:space:]]+", "%20", x),  
    '?format=geojson&limit=5&namedetails=1'),
  file = file.path(cachePath, make.names(x)),
  stringsAsFactors=FALSE
  )

  x3 = list()
  for(D in 1:nrow(xDf)) {
      if(verbose) {
        message(xDf[D,'orig'])
      }

    cacheFile = xDf[D,'file']
    if(file.exists(cacheFile)) { 
      if(verbose) {
        message("found in cache ", cachePath)
      }
    } else {
      downloadFileMapmisc(url = xDf[D,'url'], destfile = xDf[D,'file'],
        quiet=!verbose)
    }
    x3[[D]] = rgdal::readOGR(
        dsn = cacheFile,
        verbose=verbose, 
        stringsAsFactors=FALSE)

    # there might be more than one result
    # get rid of ways
    isNotWay = !x3[[D]]$osm_type == 'way'
    if(any(isNotWay)) {
        x3[[D]] = x3[[D]][isNotWay,]
      }


    # keep only places (not river, etc) if there are places
      isCat = x3[[D]]$category =='place'
      if(any(isCat)) {
        x3[[D]] = x3[[D]][isCat,]
      }
    # keep the place with the most name information
    # it's probably most important
    mostNames = which.max(nchar(x3[[D]]$namedetails))
    if(length(mostNames)) {
      x3[[D]] = x3[[D]][mostNames, ]
    }

    x3[[D]] = x3[[D]][1,]   

    named = trimws(scan(text=gsub("^[{]|[}]$", "", x3[[D]]$namedetails), 
        sep=',', what='a', quiet=TRUE))

    x3[[D]]$name = gsub(langString, "", 
      c(grep(langString, named,
      value=TRUE), NA)[1])

  }

  allNames = unique(unlist(lapply(x3, names)))
  for(D in 1:length(x3)) {
    missingHere = setdiff(allNames, names(x3[[D]]))
    for(D2 in missingHere) x3[[D]][[D2]] = NA
  }

x4=do.call(rbind, x3)
x4$orig = xDf$orig

x4$name[is.na(x4$name)] = gsub(", .*", "", x4$display_name[is.na(x4$name)])

firstCols = c('name','orig','type','category','importance')
omitCols = c('namedetails','icon')

x4@data = cbind(
  x4@data[,intersect(firstCols, names(x4))],
  x4@data[,setdiff(names(x4), c(omitCols, firstCols))]
  )

x4

}

geocodeOld = function(x, oneRecord=FALSE, extent=NULL, progress='', ...) {
  
  if(length(getOption('mapmiscVerbose'))) { 
    verbose = getOption('mapmiscVerbose')
  } else {
    verbose=FALSE
  }
  cachePath=getOption('mapmiscCachePath')
  if(is.null(cachePath)) {
    cachePath = tempdir()
  }
  if(!nchar(cachePath)) {
    cachePath = '.'
  }
  cachePath = file.path(cachePath,'geocode')
  dirCreateMapmisc(cachePath,recursive=TRUE,showWarnings=FALSE)
  
  theproj = crs(extent)
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
  }	else if(FALSE) { #requireNamespace("dismo", quietly = TRUE)) {
    if(verbose)
      message("not found in cache, retrieving")
  	Sys.sleep(1)
    result = mapmiscGeocode(
      x, oneRecord=oneRecord, 
      extent=extLL, progress=progress, ...)		
    if(verbose)
      message("caching in ", cachePath)
    
    saveRes = try(save(result, file=cacheFile), silent=TRUE)
    if(verbose & any(class(saveRes) == 'try-error')) {
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
