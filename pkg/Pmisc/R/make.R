#' Run make
#' 
#' @description Runs make with the Makefile supplied by Pmisc
#'
#' @param x a string specifying the target
#' @param suffix string to replace file extension
#' @param beamer use \code{--to=beamer}
#' @param run execute the command
#' @param ... additional arguments passed to make
#' @export
make = function(x, 
		suffix=NULL, beamer=FALSE,
		run=FALSE, ...) {

	if(!is.null(suffix))
		x = paste(
				tools::file_path_sans_ext(basename(x)),
				suffix, sep='')
	
theString= paste(' -f', 
		system.file('src','knitrMakefile', package='Pmisc'),
		x)

if(beamer)
	theString = paste(theString, 'pandocTo=beamer')

dots = list(...)
for(D in seq(from=1, by=1, len=length(dots)))
	theString = paste(theString, D)

if(run) {
	system(paste('make ', theString))
} else {
	cat(theString)
}

invisible(noquote(theString))
}
