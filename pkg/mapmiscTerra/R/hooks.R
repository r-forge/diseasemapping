
getCacheDir = function(persistent = FALSE) {
  if(persistent) {
    result = file.path(dirname(tempdir()), 
        paste("mapmiscCache", Sys.info()['user'],  sep='_')
    )
  } else {
    result = file.path(tempdir(), 'mapmiscCache')
  }
  result
}

.onAttach = function(libname, pkgname) {
  
  
  # if the cache option isn't set
  if(is.null(getOption("mapmiscCache"))) {
    
    # temporary directory
    cachePath = getCacheDir(FALSE)
    dir.create(cachePath, showWarnings=FALSE)
    options(mapmiscCachePath = cachePath)
    
    
    # check for a cache directory /tmp/mapmiscCache_user
    cachePath = getCacheDir(TRUE)
    # if it exists, set it as cache
    if(dir.exists(cachePath)) {
      if(file.access(cachePath,2)==0) {
        options(mapmiscCachePath = cachePath)
      }
    }
  }
  
  packageStartupMessage(format(paste(
              "map images will be cached in ", 
              getOption("mapmiscCachePath")))
  )
  
  invisible()
}