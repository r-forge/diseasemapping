glgm=function(data,  cells, covariates=NULL, formula=NULL, priorCI=NULL, maternRoughness=2, mesh=F,...) {
	
	# create raster for prediction
	if(!length(grep("^Raster",class(cells)))) { 
		# cells must be an integer
		cells = as.integer(cells)
		thebbox = data@bbox
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
		
	}

 	
	if(cells@nrows * cells@ncols > 10^6) warning("there are lots of cells in the prediction raster,\n this might take a very long time")
	

	# formula for inla

	names(data) = gsub("[[:punct:]]|[[:space:]]","_", names(data))
	
	if(is.null(formula))
		formula = names(data)[1]
	if(is.integer(formula))
		formula = names(data)[formula]

	if(class(formula)!= "formula") {
		if(!is.null(covariates)) {
			if(!length(names(covariates)))
				names(covariates) = paste("c", 1:length(covariates),sep="")			
			names(covariates) = gsub("[[:punct:]]|[[:space:]]","_", names(covariates))
			
		formula = as.formula(
				paste(formula, "~ ",
				paste(names(covariates),collapse="+")
			)
		)
	} else {
		formula = as.formula(paste(formula, "~1"))	
	}
	}
	allterms = rownames(attributes(terms(formula))$factors)
	notInData = allterms[! allterms %in% names(data)]
	
	if(! all(notInData %in% names(covariates)))
		warning("some terms in the model are missing from both the data and the covariates")
	for(D in notInData)
		data[[D]] = extract(covariates[[D]], data)

	
	# convert covariates to raster stack with same resulution of prediction raster.
	covariates = stackRasterList(covariates, cells)
	
	# data frame for inla
	cellsTemp = cells
	values(cellsTemp ) = NA
	data$space = extract(cellsTemp, data ,cellnumbers=T)[,"cells"]
	
	# priors
	if("sd" %in% names(priorCI)) {
		obj1 = sort(priorCI$sd^-2)
		cifun = function(pars) {
				theci = 	pgamma(obj1, shape=pars[1], rate=pars[2],log.p=T)
				
				(log(0.025) - theci[1])^2 +
				(2*(log(0.975) - theci[2]))^2		

			}
		
		precPrior2=optim(c(.5,.5/mean(obj1)), cifun, 
				lower=c(0.000001,0.0000001),method="L-BFGS-B")
		precPrior = precPrior2$par
		names(precPrior ) = c("shape","rate")
		
 		#pgamma(obj1, shape= precPrior["shape"], rate=precPrior["rate"],log.p=F)
		#pgamma(obj1, shape= precPrior["shape"], rate=precPrior["rate"],log.p=T)
		#log(c(0.025, 0.975))
		#precPrior2
		#pgamma(obj1, shape=precPrior["shape"], rate=precPrior["rate"],log.p=T)
		#log(c(0.025, 0.975)) 
 		#1/sqrt(qgamma(c(0.975,0.025), shape=precPrior["shape"], rate=precPrior["rate"]))
		#priorCI$sd
		
	} else {
		precPrior = c(shape=0.01, rate=0.01)
	}

	
	if("range" %in% names(priorCI)) {
		if(priorCI$range[1] < xres(cells)/4) {
			priorCI$range[1] = xres(cells)/4
			warning("lower bound of range CI too small, setting it to 1/4 cell size")
			
		}
		
		# rate parameter, in terms of cells, not km.
		obj1=sort(xres(cells)/priorCI$range)
		
		cifun = function(pars) {
			sum(abs(log(c(0.025, 0.975)) -
			pgamma(obj1, shape=pars[1], rate=pars[2], log.p=T)
			)*c(1, 100))
		}
			
		ratePrior2=optim(c(1,1), cifun, 
				lower=c(0.001,0.001),method="L-BFGS-B")
		ratePrior = ratePrior2$par
		names(ratePrior ) = c("shape","rate")
#		ratePrior
#		pgamma(obj1, shape=ratePrior["shape"], rate=ratePrior["rate"])
#		 qgamma(c(0.975,0.025), shape=ratePrior["shape"], rate=ratePrior["rate"])
#		xres(cells)*qgamma(c(0.975,0.025), shape=ratePrior["shape"], rate=ratePrior["rate"])
#		1/(qgamma(c(0.975,0.025), shape=ratePrior["shape"], rate=ratePrior["rate"]))
#		xres(cells)/(qgamma(c(0.975,0.025), shape=ratePrior["shape"], rate=ratePrior["rate"]))
		
	} else {
		precPrior = c(shape=0.01, rate=0.01)
	}

	spaceFormula = paste(".~.+ f(space, model='matern2d', ",
				"nrow=", nrow(cells), 
				", ncol=", ncol(cells),
				", nu=", maternRoughness, 
				", hyper = list(",
				 "range=list( param=c(",
				      paste(ratePrior, collapse=","),
				"), prior='loggamma'),",
				"prec=list( param=c(",
				paste(precPrior, collapse=","),
				"),prior='loggamma')",
				" ) )" 
			)
	
	formula = update.formula(formula,	as.formula(spaceFormula))
	
	
	# lincombs argument for prediction
	
	
	# call inla
	inlaResult = inla(formula, data=data@data, ...
	)
	#, ...)
	
	# random into raster
	
	
	# lincombs into raster
	
	
	
	# E exp(lincombs)
	
	
	# E inv logit(lincombs)
	
	
	# parameter prior and posteriors
	
		
	
	return(list(inla=inlaResult))
	
}