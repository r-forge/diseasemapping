
mapmiscCacheCheck = function() {
  if(identical(
    getOption('mapmiscCacheReadOnly'), TRUE)
  ) {
    stop("Attempting to write to read-only cache")
  }
  if(identical(
    getOption('mapmiscCachePath'), 
    system.file('extdata', package='mapmisc'))
  ) {
    stop("Attempting to write to system folder")
  }
}

downloadFileMapmisc = function(...) {
  mapmiscCacheCheck()
  utils::download.file(...)
}

dirCreateMapmisc = function(path, ...) {
  if(!dir.exists(path)) {
    mapmiscCacheCheck()
    dir.create(path, ...)
  }
}

mapmiscGeocode = function(...) {
  mapmiscCacheCheck()
#  dismo::geocode(...)
}

mapmiscGNcities = function(...) {
  mapmiscCacheCheck()
  geonames::GNcities(...)
}

persistentCache = function(verbose=TRUE) {
  
  cachePath = getCacheDir(TRUE)
  dir.create(cachePath, showWarnings=FALSE)
  
  options(mapmiscCachePath = cachePath)
  
  if(verbose)
    message(format(paste(
                "map images will be cached in ", 
                getOption("mapmiscCachePath"))
        )
    )
  
  invisible(cachePath)
  
}
