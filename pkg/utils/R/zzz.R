.onLoad = function(lib, pkg) {
	
	tools::vignetteEngine(
			'plainKnit', 
			weave=	function(file, ...) {
				# create a dummy .tex file so R CMD doesn't complain
				# make will discover that the .md file is newer and overwrite
				texFile = gsub("\\.[Rr]md$", ".tex", file)
				file.create(texFile)
				knitr::knit(input=file, ...)
			},
			tangle = knitr::purl, 
			pattern = '.[Rr]md$',
			aspell = list(
    	filter = knitr::knit_filter
			)
	)
}