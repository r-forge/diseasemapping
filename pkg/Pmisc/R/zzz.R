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
	
	if(FALSE) { # need version of R which handles #\VignetteBuilder metadata
	# spin goat's hair
	tools::vignetteEngine(
			name = 'spin', 
			package='Pmisc',
			weave =	function(file, ...) {
				knitr::spin(hair=file, format='Rmd', knit=TRUE, report=TRUE)
			},
			tangle = function(file, ...) file.copy(file, gsub("hair$", "", file)), 
			pattern = '[.][Rr]hair$'
	)
	
	# spin to an md file (don't convert to html)
	tools::vignetteEngine(
			name = 'spinReportFalse',
  		package='Pmisc',
			weave =	function(file, ...) {
				# create an empty .tex file.  
				# this is a hack because R CMD build doesn't
				# check for .md files
				file.create(gsub("\\.[[:alpha:]]+$", ".tex", file))
				knitr::spin(hair=file, format='Rmd', knit=TRUE, report=FALSE)
			},
			tangle = function(file, ...) file.copy(file, gsub("hair$", "", file)),
			pattern = '[.][Rr]hair$'
	)
}
}