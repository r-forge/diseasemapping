#' Download a file
#' 
#' @description Checks if a downloaded file is old and re-downloads is necessary
#'
#' @param url a string specifying the remote location
#' @param file local file name
#' @param age maximum age of the local file
#' @param verbose print additional information
#' @param source source the file after download
#' @export
downloadIfOld = function(
		url,
		file = basename(url),
		age = '5 days',
		verbose=FALSE, 
		source=FALSE
) {
	
	old = Sys.time() - diff(as.numeric(
					seq(from = Sys.time(), len=2, by=age)
			)) 		
	
	for(D in 1:length(url)) {
		if(verbose)
			message('file ', file[D])
		
		
		fileIsNew = any(as.numeric(difftime(
								file.info(file[D])$mtime,
								old
						)) > 0, na.rm=TRUE)
		
		if(!fileIsNew) {
			utils::download.file(
					url[D],
					file[D],
					quiet = !verbose
			)
		}
		if(source) {
			source(
					file[D],
					verbose=verbose>1
			)
		}
	}
	invisible()
}