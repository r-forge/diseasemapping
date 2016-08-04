.onLoad = function(lib, pkg) {
	
	tools::vignetteEngine(
			'plainKnit', 
			weave=	function(file, ...) {
				# create a dummy .tex file so R CMD doesn't complain
				# make will find the .md file is newer
				texFile = gsub("\\.[Rr]md$", ".tex", file)
#				pdfFile = gsub("\\.[Rr]md$", ".pdf", input)
				file.create(texFile)
				#				cat('
#\\documentclass{article}
#\\begin{document}
#stuff
#\\end{document', 
#		file=texFile)
#				outfile = knitr::knit(input, ...)
#				system(paste(
#								'pandoc --standalone --smart --to=beamer --output=',
#								texFile, input))
				knitr::knit(input=file, ...)
			},
			tangle = knitr::purl, 
			pattern = '.[Rr]md$',
			aspell = list(
    	filter = knitr::knit_filter
			)
	)
}