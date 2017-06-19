#' Convert by to data.frame
#' 
#' @description Convert the output from the by function to a data frame
#'
#' @param x produced by stats::by
#' @param format type of output, or a formula passed to dcast
#' @param name assigned to the output variable
#' @details The stats::by function produces list-like objects, which can be converted to arrays or data frames
#' @examples
#' data('meuse', package='sp')
#' byRes = by(meuse[,c('copper','lead')], meuse[, c('landuse','soil')], range)
#' 
#' dimnames(by2df(byRes, 'array'))
#' head(by2df(byRes, 'wide'))
#' head(by2df(byRes, 'long'))
#' 
#' if(requireNamespace('reshape2', quietly=TRUE)) {
#' 	head(reshape2::dcast(by2df(byRes, 'long'), landuse ~ soil + output))
#' }
#' 
#' 
#' @export
by2df = function(x, format = c('wide','long', 'array'), name='output') {
  
  theLengths = unlist(lapply(x, length))
  oneLength = setdiff(unique(theLengths), 0)
  if(length(oneLength)!= 1) warning("elements of x have different lengths")
  baseElement = min(which(theLengths != 0))  
  
  xList = lapply(x, function(xx) c(xx, rep(NA, oneLength - length(xx))))
  
  theDimnames = c(
      list(names(xList[[baseElement]])),
      dimnames(x)
      )
   names(theDimnames)[1] = name
      
   if(!length(theDimnames[[name]])) 
     theDimnames[[name]] = paste('out',1:length(xList[[baseElement]]),sep='')   
  
  res = array(unlist(xList), 
      dim = unlist(lapply(theDimnames, length)),
      dimnames = theDimnames
      )
      res = aperm(res, c(2:length(dim(res)),1))
  
  if(class(format)=='formula') {
    resDf = do.call(expand.grid, dimnames(res))
    resDf$value = as.vector(res)
    return(reshape2::dcast(resDf, formula=format))    
  }    
      
  if(format[1]=='wide') {
    resMat = matrix(res, ncol=dim(res)[length(dim(res))],
        dimnames = list(NULL, dimnames(res)[[length(dim(res))]]))
    resDf = do.call(expand.grid, dimnames(res)[seq(1,length(dim(res))-1)])
    res = cbind(resDf, resMat)    
  }      

  if(format[1]=='long') {
    resDf = do.call(expand.grid, dimnames(res))
    resDf$value = as.vector(res)
    res = resDf
  }      
    
  res
}



