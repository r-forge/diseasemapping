#' Output width for knitr figures
#' 
#' @description Creates output width strings for html and latex
#' @param x width, fraction of page width, a number between zero and 1
#' @param ... other arguments, currently ignored
#' @examples 
#' Pmisc::out.width(0.5)
#' 
#' @export
out.width = function(x, ...) {

	testPlotOutput = knitr::knit_hooks$get()$plot(
		'xx', list(fig.align='default',fig.show=TRUE))

	if(grepl("includegraphics[{]xx[}]", testPlotOutput)[1]) {
		result = paste0(x, '\\textwidth')
	} else {
		result = paste0(100*x, '%')
	}
	
	result
}
