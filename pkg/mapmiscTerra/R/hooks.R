
getCacheDir = function(persistent = FALSE) {
  # Under R CMD check, always use the session tempdir so a persistent
  # cache next to it is not left as temp detritus.
  inCheck = nzchar(Sys.getenv("_R_CHECK_PACKAGE_NAME_", unset = ""))
  if(persistent && !inCheck) {
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
    
    
    # Prefer an existing persistent cache when available (no-op under R CMD check).
    cachePath = getCacheDir(TRUE)
    if(dir.exists(cachePath) && !identical(cachePath, getOption("mapmiscCachePath"))) {
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