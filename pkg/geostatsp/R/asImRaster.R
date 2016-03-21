#as.im = function(X, ...) {
#	UseMethod("as.im")
#}

asImRaster = function(X, ...) {
  if(requireNamespace('spatstat',quietly=TRUE)){	
	
	if(nlayers(X)>1) X = X[[1]]	
	xMat = raster::as.matrix(X)[nrow(X):1,]
	if(is.factor(X)) {
		facTable = levels(X)[[1]]
		xMat = matrix(
				facTable[match(xMat, facTable[,1]), min(c(2,ncol(facTable)))], 
				nrow(xMat), ncol(xMat))
		levels(xMat) = facTable[,min(2, ncol(facTable))]
	}
		
	res=spatstat::as.im(xMat, 
			xrange=bbox(X)[1,], 
			yrange=bbox(X)[2,], ...)
	
} else {
	message("Install the spatstat package if you wish to use this function.")
	res = NULL
}
return(res)
}

#as.im.default = function(X, ...) {
#	spatstat::as.im(X, ...) 
#}