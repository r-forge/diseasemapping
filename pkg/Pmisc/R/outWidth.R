#' Output width for knitr figures
#' 
#' @description Creates output width strings for html and latex
#' @param x width, fraction of page width
#' @param mdToTex if \code{FALSE} convert to markdown table, return \code{x} otherwise
#' @examples 
#' Pmisc::out.width(0.5, 'latex')
#' Pmisc::out.width(0.5, 'markdown')
#' mdToTex=FALSE
#' Pmisc::out.width(0.5, 'auto')
#' mdToTex=TRUE
#' Pmisc::out.width(0.5, 'auto')
#' 
#' @export
out.width = function(x, mdToTex = 'auto') {
	if(identical(mdToTex, 'auto')) {
		mdToTex = any(commandArgs()=='mdToTex', na.rm=TRUE)
		if(exists('mdToTex')) {
			mdToTex = identical(get("mdToTex"), TRUE)
		}
	}
	
	if(identical(mdToTex, TRUE) | identical(mdToTex, 'latex') ) {
		
		result = paste(x, '\\textwidth', sep='')
	} else {
		result = paste(100*x, '%', sep='')
	}
	
	result
}
