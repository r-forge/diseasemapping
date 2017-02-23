.onLoad = function(lib, pkg) {
	
	# basic knitr builder
	tools::vignetteEngine(
			name='barebones', 
			package='Pmisc',
			weave= function(file, ...) {
				# create an empty .tex file.  
				# this is a hack because R CMD build doesn't
				# check for .md files
				file.create(gsub("\\.[[:alpha:]]+$", ".tex", file))
				knitr::knit(file, ...)
			},
			tangle = function(...) knitr::purl(...), 
			pattern = "[.][rR](nw|md)$",
			aspell = list(
    			filter = function(...) knitr::knit_filter(...)
			)
	)
	
	# spin goat's hair
	tools::vignetteEngine(
			name = 'spin', 
			package='Pmisc',
			weave =	function(file, ...) {
				knitr::spin(hair=file, format='Rmd', knit=TRUE, report=TRUE)
			},
			tangle = function(file, ...) {
				hair=readLines(file)
				wool=knitr::spin(
						text=hair, 
						format='Rmd', knit=FALSE)
				knitr::knit(
						text=wool,
						tangle=TRUE, 
						output = gsub("hair$", "", file)	
				)
			}, 
			pattern = '[.][Rr]hair$'
	)
	
}
