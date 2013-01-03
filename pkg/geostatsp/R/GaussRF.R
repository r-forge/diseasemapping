GaussRF = function(x, ...) {
	UseMethod("GaussRF")
	
}

GaussRF.Raster = function(x, ...){

	res =RandomFields::GaussRF(
			x = seq(x@extent@xmin, x@extent@xmax, len=x@ncols),
			y = seq(x@extent@ymin, x@extent@ymax, len=x@nrows),
			grid=TRUE, 
		...
			)

			
	resRast = raster(t(res[,seq(dim(res)[2], 1)]),
			x@extent@xmin, x@extent@xmax,
			x@extent@ymin, x@extent@ymax, crs=x@crs)
	
	return(resRast)
}

GaussRF.SpatialPointsDataFrame = function(x, ...){
	
	x$GaussRF =RandomFields::GaussRF(
		x@coords,
		...
	)

	return(x)
}

GaussRF.SpatialPoints= function(x, ...){
	
	x = SpatialPointsDataFrame(x, data=data.frame(GaussRF= rep(NA, length(x))))
	x$GaussRF =RandomFields::GaussRF(
			x@coords,
			...
	)
	
	return(x)
}

GaussRF.default = function(x,  ...){
	

	RandomFields::GaussRF(x,
			...
	)
	
	
}


