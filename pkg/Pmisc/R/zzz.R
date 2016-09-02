.onLoad = function(lib, pkg) {
	
	# basic knitr builder
	tools::vignetteEngine(
			name='barebones', 
			package='Pmisc',
			weave= function(...) knitr::knit(...),
			tangle = function(...) knitr::purl(...), 
			pattern = c(
					input="[.]([rRsS](nw|tex)|[Rr](md|html|rst))$", 
					output="[.](tex|md|html|pdf|rst)$"),
			aspell = list(
    			filter = function(...) knitr::knit_filter(...)
			)
	)
	
	# spin goat's hair
	tools::vignetteEngine(
			name = 'spin', 
			package='Pmisc',
			weave =	function(file, ...) knitr::spin(file),
			tangle = function(file, ...) file.copy(file, gsub("hair$", "", file)), 
			pattern = '[.][Rr]hair$'
	)
	
	# spin without converting md to html
	tools::vignetteEngine(
			name = 'spinReportFalse',
  		package='Pmisc',
			weave =	function(file, ...) {
				knitr::spin(
						hair=file, 
						format='Rmd', knit=TRUE, 
						report=FALSE, precious=TRUE)
			},
			tangle = function(file, ...) file.copy(file, gsub("hair$", "", file)), 
			pattern = c(input='[.][Rr]hair$',output='[.]md$')
	)
	
}