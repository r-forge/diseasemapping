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
	
	result =RandomFields::GaussRF(
		x@coords,
		...
	)

	return(result)
}

GaussRF.SpatialPoints= function(x, ...){
	
	result =RandomFields::GaussRF(
			x@coords,
			...
	)
	
	return(result)
}

GaussRF.default = function(x,  ...){
	

	RandomFields::GaussRF(x,
			...
	)
	
	
}


