lgcp = function(data,  cells, covariates=NULL, formula=NULL, priorCI=NULL, 
rough=1, buffer = 0, mesh=F,...) {

# create raster for prediction
if(!length(grep("^Raster",class(cells)))) { 
	# cells must be an integer
	cells = as.integer(cells)
	smallBbox = thebbox = data@bbox
	thebbox = thebbox + buffer*cbind(-c(1,1),c(1,1))
	res = diff(thebbox[1,])/cells		
	Nx = cells
	Ny = ceiling(diff(thebbox[2,])/res)
	thebbox[2,2] = thebbox[2,1] + res*Ny
	cells= raster(extent(thebbox), ncols=Nx, nrows=Ny, crs=data@proj4string)
} else {
	# it's a raster, make sure it has square cells
	if(xres(cells) != yres(cells)) 
		res = xres(cells)
	theextent = cells@extent
	theylim = theextent@ymax - theextent@ymin
	Ny = ceiling(theylim/res)
	theextent@ymax = theextent@ymin + Ny * res
	
	cells = raster(theextent, ncols=cells@ncols, nrows=Ny,crs=cells@crs)
	smallBbox = bbox(cells)
	
}

# create data


	
	data = rasterize(data, cells, fun="count")
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
	formula.tools::lhs(formula) = as.name("count")

	
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
		formula=formula,priorCI=priorCI,rough=rough,
		buffer=buffer, mesh=mesh, 
		family="poisson",
		...)

result

}

if(F) {

	mesh=F; 
	family="poisson"
	
}

