simPoissonPP = function(intensity) {
	
	intensity = intensity * prod(res(intensity))
	
	tosub = data.frame(id=NA,v=0)
	intensity = subst(intensity, tosub$id, tosub$v)

  if(any(minmax(intensity)>1000))
    warning("A large number of events are being simulated, more than ", max(minmax(intensity)))
  
  
	NperCell = intensity
	terra::values(NperCell)= rpois(
      ncell(intensity)*nlyr(intensity), 
      values(intensity))
	
	if(any(minmax(NperCell)>1000))
		warning("A large number of events are being simulated, more than", max(minmax(NperCell)))
	
	events= vector('list', nlyr(intensity))
	for(D in 1:nlyr(intensity)) {
	eventsD = rep(1:ncell(NperCell), values(NperCell[[D]]))
	
	if(length(eventsD)>1e6)
		warning("more than 1,000,000 events being simulated")
	
	eventsD = xyFromCell(NperCell, eventsD) 
	
	eventsD= eventsD + cbind(
			runif(dim(eventsD)[1],-xres(intensity)/2, xres(intensity)/2),
			runif(dim(eventsD)[1],-yres(intensity)/2, yres(intensity)/2)
	)
	
  if(nrow(eventsD)){
  	events[[D]] = vect(eventsD)
	  terra::crs(events[[D]]) = crs(intensity)	
  } else {
    events[[D]] = NULL
  }
	} # D loop

	if(length(events)==1) {
		names(events) = 'events'
	} else {
		names(events) = paste("events", 1:length(events),sep="")	
	}
	events
}

simLgcp = function(param, covariates=NULL, betas=NULL, 
		offset=NULL, 
		rasterTemplate=covariates[[1]],  n=1, 
		...) {
	
	randomEffect = RFsimulate(model=param, x=rasterTemplate, n=n, ...)

	if(length(offset)==length(names(covariates)))
    offset = names(covariates)
	
	if(is.null(names(betas)) & (!length(offset)))
		names(betas) = names(covariates)
	
	if(!is.null(covariates)){
		# method = bilinear?
		covariates = stackRasterList(covariates, randomEffect)
	}
	
	if(length(betas) & is.null(names(betas))) {
		names(betas) = names(covariates)[1:length(betas)]
	}
	
	betas = c(rep(1, length(offset)), betas)
	names(betas)[seq(1, len=length(offset))] = offset

		
	covariates = covariates[[intersect(names(covariates), 
					names(betas))]]
	
	themean = 0
	if('mean' %in% names(param))
		themean = themean + param['mean']
	if('intercept' %in% names(betas))
		themean = themean + betas['intercept']
	betas['intercept'] = themean
	param = param[! names(param) %in% "mean"]

	
	
	thefixed = rast(randomEffect)
  terra::values(thefixed) = themean
  for(Dbeta in names(covariates))
    thefixed = thefixed + betas[Dbeta]*covariates[[Dbeta]]

  linearPredictor = thefixed[[rep(1, n)]]
	linearPredictor = linearPredictor +randomEffect 
	names(linearPredictor) = gsub("^sim", "linearPredictor", names(randomEffect))
	
	intensity = exp(linearPredictor)
	names(intensity) = gsub("^sim", "intensity", names(randomEffect))


	intSansOffset = linearPredictor
	for(Doffset in offset){
		intSansOffset = intSansOffset - covariates[[Doffset]]
	}
	intSansOffset = exp(intSansOffset)
	names(intSansOffset) = gsub("^sim", "relativeIntensity", names(randomEffect))
	
	events = simPoissonPP(intensity)
	
	names(randomEffect) = gsub("^sim", "random", names(randomEffect))
	
	return(c(
			events, 
			list(
				raster = c(randomEffect,
						linearPredictor, intensity,
						intSansOffset,
						covariates),
				parameters=list(random=param, fixed=betas)
			)
		)
	)
}