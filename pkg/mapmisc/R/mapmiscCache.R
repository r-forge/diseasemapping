


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