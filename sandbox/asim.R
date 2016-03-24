as.im.RasterLayer = function (from) {
  if (!requireNamespace("spatstat", quietly = TRUE)) 
    stop("package spatstat required for coercion")
  if (!requireNamespace("raster", quietly = TRUE)) 
    stop("package raster required for coercion")
  if (!raster::hasValues(from)) 
    stop("values required in RasterLayer object")
  if (raster::rotated(from)) {
    stop("\n Cannot coerce because the object is rotated.\n Either coerce to SpatialPoints* from\n or first use the \"rectify\" function")
  }
	
	if(is.factor(from)) {
		# create an im with levels
		facTable = levels(from)[[1]]
		xMat = raster::as.matrix(from)[nrow(from):1,]
		
		xMat = matrix(
				facTable[match(xMat, facTable[,1]), min(c(2,ncol(facTable)))], 
				nrow(xMat), ncol(xMat))
		
		levels(xMat) = facTable[,min(2, ncol(facTable))]
		
		im=spatstat::as.im(xMat, 
			xrange=bbox(from)[1,], 
			yrange=bbox(from)[2,])
	
	} else {
	
  	rs <- raster::res(from)
  	orig <- bbox(from)[, 1] + 0.5 * rs
  	dm <- dim(from)[2:1]
  	xx <- unname(orig[1] + cumsum(c(0, rep(rs[1], dm[1] - 1))))
  	yy <- unname(orig[2] + cumsum(c(0, rep(rs[2], dm[2] - 1))))
  	im <- spatstat::im(matrix(raster::values(from), ncol = dm[1], 
        	nrow = dm[2], byrow = TRUE)[dm[2]:1, ], xcol = xx, yrow = yy)
	}
  im
}

example = function() {
	
	# numeric raster
	myRaster = raster(matrix(100+1:5,10,10),
			xmn=1,xmx=6,ymn=-1,ymx=4)
	
	myIm = as.im.RasterLayer(myRaster)

	# factor raster
	myRasterF = myRaster
	levels(myRasterF) = data.frame(
			ID = sort(unique(values(myRaster))), 
			labels=letters[seq(1, length(unique(values(myRaster))))]
	)
	myImF = as.im.RasterLayer(myRasterF)

	myIm
	myImF

}