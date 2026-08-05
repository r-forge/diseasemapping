
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

downloadFileMapmisc = function(url, destfile, ...) {
  mapmiscCacheCheck()
  dots = list(...)
  # OSM tile usage policy requires an identifying User-Agent
  ua = paste0(
    "mapmisc/",
    tryCatch(as.character(utils::packageVersion("mapmisc")),
      error = function(e) "devel"),
    " (https://diseasemapping.r-forge.r-project.org/; R package)"
  )
  headers = dots$headers
  if(is.null(headers)) headers = character()
  if(!length(grep("(?i)^user-agent$", names(headers)))) {
    headers = c(headers, "User-Agent" = ua)
  }
  dots$headers = headers
  dots$url = url
  dots$destfile = destfile
  do.call(utils::download.file, dots)
}

dirCreateMapmisc = function(path, ...) {
  if(!dir.exists(path)) {
    mapmiscCacheCheck()
    dir.create(path, ...)
  }
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
