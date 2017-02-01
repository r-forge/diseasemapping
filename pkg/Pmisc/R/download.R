#' Download a file
#' 
#' @description Checks if a downloaded file is old and re-downloads is necessary
#'
#' @param url a string specifying the remote location
#' @param file local file name
#' @param path local directory to store downloaded files
#' @param age maximum age of the local file
#' @param verbose print additional information
#' @param ... additional arguments for \code{\link[utils]{download.file}}
#' @return A character string of files downloaded or unzipped.
#' @details Checks if a downloaded file is old and re-downloads is necessary.  Any zip files are unzipped.
#' @examples
#' \dontrun{  
#'   theFiles = Pmisc::downloadIfOld(
#'     'https://cran.r-project.org/src/base/Historic/Windows/mva.zip'
#'   )
#'   theFiles
#' }

#' @seealso \code{\link[base]{download.file}}, \code{\link[utils]{unzip}}
#' @export
downloadIfOld = function(
  url,
  path = getwd(),
  file = file.path(path, basename(url)),
  age = '5 days',
  verbose=FALSE,
  ...
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
        quiet = !verbose, 
        ...
      )
    }
    
  }
  
  
  gzFiles = grep("[.]gz$", file, ignore.case=TRUE)  
  if(requireNamespace('R.utils')) {
    for(D in gzFiles) {
      outfile = gsub("[.]gz$", "", file[D], ignore.case=TRUE)
      fileIsNew = any(as.numeric(difftime(
            file.info(outfile)$mtime,
            old
          )) > 0, na.rm=TRUE)
      if(!fileIsNew) {        
        file[D] = R.utils::gunzip(
          file[D], 
          overwrite = file.exists(gsub(
              "[.]gz$", "", file[D], ignore.case=TRUE
            )),
          remove=FALSE
        )
      } else {
        file[D]= outfile
      }
    }
  }
  
  res = file
  exdir = rep_len(path, length(file))
  
  zipfiles = grep('[.]zip$', file, ignore.case=TRUE)
  
  for(D in zipfiles) {
    res = c(res, utils::unzip(file[D], exdir=exdir[D]))
  }
  
  invisible(res)
  
}