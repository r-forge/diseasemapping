lgcp = function(data,  cells, covariates=NULL, 
		formula=NULL, priorCI=NULL, 
shape=1, buffer = 0, mesh=FALSE,...) {

	
	# create raster for prediction
	if(!length(grep("^Raster",class(cells)))) { 
		# cells must be an integer
		cells = squareRaster(data, cells)
	} else {
		cells = squareRaster(cells)
	}
	smallBbox = bbox(cells)
	buffer =  ceiling(buffer/xres(cells))
	cells = extend(cells, c(buffer, buffer))
	thebbox = bbox(cells)	
	

# create data


	
	data = rasterize(SpatialPoints(data), cells, fun="count")
	names(data) = "count"
	data[is.na(data)] = 0
	
# the formula	
	if(is.null(formula)) {
		formula = as.formula(
				paste(c("count ~ 1", names(covariates)), collapse="+")
		)
	}

	formula	= update.formula(formula,
			.~.+offset(logCellSize) 
	)
	formula = update.formula(formula, count ~ .)

	
	# cell size offset
	logCellSize = cells
	names(logCellSize) = "logCellSize"
	values(logCellSize) =  sum(log(res(cells)) )

  	if(class(covariates)=="RasterLayer") {
		covariates = list( logCellSize, covariates)
		names(covariates) = unlist(lapply(covariates, names))
	} else {
		if(is.null(covariates)) covariates = list()
		covariates = c(covariates, logCellSize=logCellSize)
	}

 result = glgm(data=data, cells=cells, covariates=covariates, 
		formula=formula,priorCI=priorCI,shape=shape,
		buffer=buffer, mesh=mesh, 
		family="poisson",
		...)

result

}


