.onLoad = function(lib, pkg) {
	
	# basic knitr builder
	tools::vignetteEngine(
			name='barebones', 
			package='Pmisc',
			weave= function(...) knitr::knit(...),
			tangle = function(...) knitr::purl(...), 
			pattern = "[.]([rRsS](nw|tex)|[Rr](md|html|rst))$",
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
			tangle = function(file, ...) file.copy(file, gsub("hair$", "", file)), 
			pattern = '[.][Rr]hair$'
	)
	
	# spin to an md file (don't convert to html)
	tools::vignetteEngine(
			name = 'spinReportFalse',
  		package='Pmisc',
			weave =	function(file, ...) {
				knitr::spin(hair=file, format='Rmd', knit=TRUE, report=FALSE)
			},
			tangle = function(file, ...) file.copy(file, gsub("hair$", "", file)),
			pattern = '[.][Rr]hair$'
	)
	
}