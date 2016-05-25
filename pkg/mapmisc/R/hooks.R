.onLoad = function(lib, pkg) {
	# check for cachePath
	if(!is.null(options()$mapmisc['cachePath'])) {
		packageStartupMessage("map images will be cached in ", options()$mapmisc['cachePath'])		
	}
	if(!is.null(options()$mapmisc['verbose'])) {
		packageStartupMessage("verbose is ", options()$mapmisc['verbose'])		
	}
	invisible()
}
