simLgcp = function(param, covariates=NULL, betas=NULL, 
		rasterTemplate=covariates[[1]], model="whittle", ...) {
	
	param["nugget"] = 0
	param = param[c("mean","variance","nugget", 
				"scale" ,"alpha")]
	
	if(!is.null(covariates))
		covariates = stackRasterList(covariates, rasterTemplate)

	randomEffect = GaussRF(rasterTemplate, model=model, 
					param=param, ...)

	linearPredictor = randomEffect
	if(is.null(names(betas)))
		names(betas) = names(covariates)
	
	for(D in names(covariates)) {
		linearPredictor = linearPredictor + betas[D]* covariates[[D]]		
	}

	intensity = exp(linearPredictor)
	
	NperCell = intensity
	values(NperCell) = rpois(ncell(intensity), 
			values(intensity)*prod(res(rasterTemplate)))
	
	events = rep(1:ncell(NperCell), values(NperCell))
	
	events = as.data.frame(NperCell,xy=TRUE)[events, c("x","y")]
	
	events = events + cbind(
			runif(dim(events)[1],-xres(rasterTemplate)/2, xres(rasterTemplate)/2),
			runif(dim(events)[1],-yres(rasterTemplate)/2, yres(rasterTemplate)/2)
	)
	
	events = SpatialPoints(events)
	
	if(!is.na(rasterTemplate@crs@projargs))
		events@proj4string = rasterTemplate@crs
	
	names(linearPredictor) = "linearPredictor"
	names(intensity) = "intensity"
	names(randomEffect) = "random"
	
	return(list(events=events,
					raster = stack(randomEffect,
							linearPredictor, intensity,covariates),
					params=list(random=param, fixed=betas)
			))
			

}