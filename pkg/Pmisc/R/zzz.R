.onLoad = function(lib, pkg) {
	
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
