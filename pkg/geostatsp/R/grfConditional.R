grfConditional = function(data, y=1, 
			param, locations, Nsim,
		 	fun=NULL, nuggetInPrediction=TRUE){
		
		
	if(is.numeric(locations)){
		# locations is number of cells in the x direction
		Nx = locations[1]
		myExtent = 	extent(data@bbox)
		Ny = round(locations*diff(data@bbox[2,])/diff(data@bbox[1,]))
		myExtent@ymax = myExtent@ymin + Ny * diff(data@bbox[1,])/Nx
		locations = raster(myExtent, Ny, Nx,
				crs=projection(data))	
	}
	if(nrow(locations) * ncol(locations) > 10^7) warning("there are lots of cells in the prediction raster,\n this might take a very long time")
	
	xseq = c(xmin(locations)+xres(locations)/2, 
			xmax(locations)-xres(locations)/2, xres(locations))
	yseq = c(ymin(locations)+yres(locations)/2,
			ymax(locations)-yres(locations)/2, yres(locations))


	if(length(y)==1 & length(names(data)) > 1) {
		y = data.frame(data)[,y]
	}

	
	if(length(dim(y))==3) {
		y = matrix(y, nrow=prod(dim(y)[1:2]), ncol=dim(y)[3])
	}

	# convert param to a matrix if it's a list
	# fill in missing parameters with defaults
	param = fillParam(param)
	if(class(param)%in%c("integer","numeric"))
		param = t(as.matrix(param))
	if(class(y)%in%c("integer","numeric")) {
		y = t(as.matrix(y))
	
	Nsamples = unique(c(dim(param)[1], dim(y)[1]))
	Nsamples = Nsamples[Nsamples!=1]
	if(length(Nsamples)>1.5)
		warning("number of samples in y and param is different")
	if(!length(Nsamples))
		Nsamples = 1
	
	param = param[round(seq(1,dim(param)[1], len=Nsim)),,drop=FALSE]
	y = y[round(seq(1,dim(y)[1], len=Nsim)),,drop=FALSE]
}	




resTemplate = raster(locations)

simFun = function(D) {

	if(param[D,"nugget"] > 0) {
		err.model = "nugget"
		err.param=c(1,param[D,"nugget"],0)
	} else {
		err.model = err.param = NULL
	}
	
	modelv = modelRandomFields(param[D,])
	
	res = RandomFields::CondSimu(krige.method="O",
			x=xseq, y = yseq, grid=TRUE, gridtriple=TRUE,
			param=NULL, model=modelv,
			given=data@coords,
			data=y[D,],
			err.model=err.model,
			err.param=err.param, method="direct decomp."
	)		
	
	
#	CondSimu("S", given=locations.obs, 
#			data=params[[theEffectR]][Siter[Diter],Dchain,], 
#			x=xgrid, y=ygrid, grid=TRUE, model="exponential", 
#			param=c(mean=0, 
#					variance=params[[theSD]][Siter[Diter],Dchain]^2, 
#					nugget=0, 
#					scale=params[[thePhi]][Siter[Diter],Dchain]), 
#			pch=" ")
	
	values(resTemplate) = t(res[nrow(res):1,])
	res = resTemplate	
	
	if(nuggetInPrediction){
		values(resTemplate) = 
				rnorm(ncell(res), sd=sqrt(param[,"nugget"]))
		res = res + resTemplate
	}		
	if(!is.null(fun)) {
		res = fun(res)
	}

	res
	}		

	result = mcmapply(simFun, 1:Nsim, SIMPLIFY=TRUE)
	if(all(sapply(result, class)=="RasterLayer"))
		result = do.call(brick, result)

result	
}


