.onAttach = function(libname, pkgname) {
	# check for cachePath

	if(.Platform$OS.type == 'unix' & is.null(options()$mapmiscCachePath)) {
		# create a cache director /tmp/user_mapmiscCache
		
		cachePath = file.path(
				gsub("/Rtmp[[:alnum:]]+$", "", tempdir()),
				paste("mapmiscCache", Sys.info()['user'],  sep='_')
		)
		dir.create(cachePath, showWarnings=FALSE)
		if(file.access(cachePath,2)==0) {
			options(mapmiscCachePath = cachePath)
		}
	}
	
	if(!is.null(options()$mapmiscCachePath)) {
		packageStartupMessage(format(paste(
			"map images will be cached in ", options()$mapmiscCachePath))
		)
	} 
	invisible()
}
