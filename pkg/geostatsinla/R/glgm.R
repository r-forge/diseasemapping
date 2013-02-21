glgm=function(data,  cells, covariates=NULL, formula=NULL, 
		priorCI=NULL, maternRoughness=2, buffer=0,
		mesh=F,...) {
	
	# create raster for prediction
	if(!length(grep("^Raster",class(cells)))) { 
		# cells must be an integer
		cells = as.integer(cells)
		thebbox = data@bbox
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
			theci = 		pgamma(obj1, shape=pars[1], rate=pars[2], log.p=T)
			
			(theci[1] - log(0.025))^2 +
					(theci[2] - log(0.925))^2 
		}
			
		ratePrior2=optim(c(2,2/mean(obj1)), cifun, 
				lower=c(0.001,0.001),method="L-BFGS-B")
		ratePrior = ratePrior2$par
		names(ratePrior ) = c("shape","rate")
if(F) {
		ratePrior
		pgamma(obj1, shape=ratePrior["shape"], rate=ratePrior["rate"])
		 qgamma(c(0.975,0.025), shape=ratePrior["shape"], rate=ratePrior["rate"])
		xres(cells)/qgamma(c(0.975,0.025), shape=ratePrior["shape"], rate=ratePrior["rate"])
	}
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
	thelincombs = as.data.frame(covariates)
	thelincombs[,"(Intercept)"]=1
	thelincombs = inla.make.lincombs(thelincombs)
	cellIndex = length(thelincombs[[1]])+1
	for(D in 1:length(thelincombs)) {
		thelincombs[[D]][[cellIndex]] = 
				list(
					space=list(idx=D,weight=1))
	}
	names(thelincombs) = paste("c", 1:length(thelincombs),sep="")
	
	# call inla
	inlaResult = inla(formula, data=data@data,
			lincomb=thelincombs,
	...
	)

	# parameter priors 
	
	params = list(
			range = list(userPriorCI=priorCI$range, 
					priorCI = 
							xres(cells)/
							qgamma(c(0.975,0.025), 
									shape=ratePrior["shape"], 
									rate=ratePrior["rate"]),
					priorCIcells = 
							1/
							qgamma(c(0.975,0.025), 
									shape=ratePrior["shape"], 
									rate=ratePrior["rate"]),
					params.intern = ratePrior),
			sd = list(userPriorCI=priorCI$sd, 
					priorCI = 
							1/sqrt(
									qgamma(c(0.975,0.025), 
											shape=precPrior["shape"], 
											rate=precPrior["rate"])),
					params.intern=precPrior)
	)
	
	rateLim = 	qgamma(c(0.999,0.001), 
			shape=ratePrior["shape"], 
			rate=ratePrior["rate"])
	rangeLim = xres(cells) /rateLim
	rangeSeq = seq(min(rangeLim), max(rangeLim), len=1000)
	rateSeq = xres(cells)/rangeSeq
	params$range$priorDist=cbind(
			x=rangeSeq,
			y=dgamma(rateSeq, shape=ratePrior["shape"], 
					rate=ratePrior["rate"]) * (rangeSeq)^(-2) * xres(cells)
		)
	if(F) {
		plot(params$range$priorDist, type='l')
	
		sum(params$range$priorDist[,"y"])*range(diff(params$range$priorDist[,"x"]))
	}	
	
	precLim = 	qgamma(c(0.999,0.001), 
			shape=precPrior["shape"], 
			rate=precPrior["rate"])
	sdLim = 1/sqrt(precLim)
	sdSeq = seq(min(sdLim), max(sdLim), len=1000)
	precSeq = sdSeq^(-2)
	params$sd$priorDist=cbind(
			x=sdSeq,
			y=dgamma(precSeq, shape=precPrior["shape"], 
					rate=precPrior["rate"]) *2* (precSeq)^(3/2) 
	)
	if(F) {
		plot(params$sd$priorDist, type='l')
		
		sum(params$sd$priorDist[,"y"])*range(diff(params$sd$priorDist[,"x"]))
	}	
	
	
	
	# random into raster
	resRasterRandom = 
			brick(cells,
			nl=dim(inlaResult$summary.random[["space"]])[2]-1)
	if("summary.random" %in% names(inlaResult)) {
		resRasterRandom[] = as.matrix(inlaResult$summary.random[["space"]][,-1])
		names(resRasterRandom) = paste("random.", names(resRasterRandom),sep="")
	} else {
		return(list(inla=inlaResult, parameters=params))
	}
	

	# E exp(lincombs)
	temp=unlist(
			lapply(inlaResult$marginals.lincomb.derived, function(qq) {
						sum(
								qq[,"x"]*c(0,diff(exp(qq[,"x"])))*qq[,"y"]	
						)
					})
	)
	inlaResult$summary.lincomb.derived[,"mean.exp"] = temp

	# E exp(lincombs)
	temp=unlist(
			lapply(inlaResult$marginals.lincomb.derived, function(qq) {
						sum(
								exp(qq[,"x"])*c(0,diff(qq[,"x"]))*qq[,"y"]	
						)
					})
	)
	inlaResult$summary.lincomb.derived[,"exp"] = temp
	
	# E inv logit(lincombs)
	if(length(grep("binomial",inlaResult$.args$family))) {
		temp=unlist(
				lapply(inlaResult$marginals.lincomb.derived, function(qq) {
							eqqx = exp(qq[,"x"])
							sum(
								eqqx/(1+eqqx)*c(0,diff(qq[,"x"]))*qq[,"y"]	
							)
						})
		)
		inlaResult$summary.lincomb.derived[,"invlogit"] = temp		
	}

		
	# lincombs into raster
	resRasterFitted = 
			brick(cells,
					nl=dim(inlaResult$summary.lincomb.derived)[2]-1)
	
		resRasterFitted[] = as.matrix(inlaResult$summary.lincomb.derived[,-1])
		names(resRasterFitted) = paste("predict.", names(resRasterFitted),sep="")
	

	# posterior distributions
parameters$sd$posterior=NULL
parameters$range$posterior=NULL

parameters$summary = inlaResult$summary.fixed
parameters$summary = rbind(parameters$summary,
		sd=NA, range=NA)
		
	
	result=list(inla=inlaResult,
					raster=stack(resRasterRandom, resRasterFitted),
					parameters=params
			)

	result
	
}