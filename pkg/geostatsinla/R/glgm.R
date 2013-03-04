glgm=function(data,  cells, covariates=NULL, formula=NULL, 
		priorCI=NULL, maternRoughness=1, buffer=0,
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
 		if(diff(res(cells))>10^(-6))  {
			res = xres(cells)
			theextent = cells@extent
			theylim = theextent@ymax - theextent@ymin
			Ny = ceiling(theylim/res)
			theextent@ymax = theextent@ymin + Ny * res
			
			cells = raster(theextent, ncols=cells@ncols, nrows=Ny,crs=cells@crs)
		}
	}
 
	
 	if(cells@nrows * cells@ncols > 10^6) warning("there are lots of cells in the prediction raster,\n this might take a very long time")
	

	
	# the formula
	# get rid of special character is names of data
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
	# convert covariates to raster stack with same resulution of prediction raster.
	if(!is.null(covariates)){
		method = rep("ngb", length(covariates))
		covariates = stackRasterList(covariates, cells, method=method)
		
	} 

	# create data frame for inla
	# if data is a raster
	if(length(grep("^Raster", class(data)))) {
 
	data = stack(data, covariates)

	data=as.data.frame(data)
		data$space = seq(1, dim(data)[1])
		
		# get rid of any rows with -Inf or NA, 
		# usually offset of zero on the natural scale
		isMinusInf = apply(data, 1, function(qq) any(qq==-Inf | is.na(qq) ))
		data = data[!isMinusInf,]
		
	}  else { # data is a SpatialPointsDataFrame
		
		notInData = allterms[! allterms %in% names(data)]
		
		if(! all(notInData %in% names(covariates)))
			warning("some terms in the model are missing from both the data and the covariates")
		for(D in notInData)
			data[[D]] = extract(covariates[[D]], data)
		
		cellsTemp = cells
		values(cellsTemp ) = NA
		data$space = extract(cellsTemp, data ,cellnumbers=T)[,"cells"]
		data = data@data
	}
	# data is now a data frame.
	
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
		ratePrior = c(shape=0.01, rate=0.01)
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
	
	
	# create linear combinations object for prediction.
	# create formula, strip out left variable and f(...) terms
	formulaForLincombs = strsplit(as.character(formula), "~")[[1]][2]
	formulaForLincombs =
			gsub("\\+?[[:space:]]*f\\([[:print:]]*\\)[[:space:]]?($|\\+)", "+", formulaForLincombs)
	# strip out offsets
formulaForLincombs =
		gsub("\\+?[[:space:]]*offset\\([[:print:]]*\\)[[:space:]]?($|\\+)", "+", formulaForLincombs)

	# convert multiple + to a single +
	formulaForLincombs = gsub(
			"\\+[[:space:]]?\\+([[:space:]]?\\+)?", "+",
			formulaForLincombs)
	# strip out trailing +
formulaForLincombs = gsub("\\+[[:space:]]?$", "", formulaForLincombs)

	if(nchar(formulaForLincombs)) {
		formulaForLincombs=as.formula(paste("~", formulaForLincombs))
		
	# model matrix from covariates	
	thelincombs = as.data.frame(covariates)
	someMissing = apply(thelincombs, 1, function(qq) any(is.na(qq)))
	# variables in the model but not in prediction rasters
	thevars = rownames(attributes(terms(formulaForLincombs))$factors)
	varsInPredict = thevars[thevars %in% names(thelincombs)]
	cantPredict = thevars[! thevars %in% names(data)]
	theFactors = grep("factor\\([[:print:]]*\\)", cantPredict)
	if(length(theFactors)) {
		temp = cantPredict
		cantPredict = cantPredict[-theFactors]
		theFactorsInFormula = temp[theFactors]
	}
	thelincombs[,cantPredict]= 0

	# check for factors	
	thefactors = unlist(lapply(data[,varsInPredict,drop=F], is.factor))
	# convert raster data into factors using same levels as in points data
	if(any(thefactors)) {
		for(D in varsInPredict[thefactors]) {
			thetable = table(data[,D])

			baseline = as.integer(names(thetable)[min(which(thetable > 0))])
			# levels not in the points data changed to baseline

			dontHave = ! (thelincombs[,D] %in% as.integer(names(thetable)))
			
			thelincombs[   dontHave   ,D] = baseline
			
			# replace NA's with baseline

			thelincombs[ is.na(thelincombs[,D] )  ,D] = baseline
			
			thelincombs[,D] = factor(thelincombs[,D],
					levels=levels(data[,D]))
		}
	}
	
	thelincombs[is.na(thelincombs)] = 0
	lincombMat = model.matrix(formulaForLincombs, thelincombs, na.action=NULL)
#	return(list(lincombMat, thelincombs, covariates, formulaForLincombs))

} else {
	lincombMat = matrix(rep(1,ncell(cells)), ncell(cells),1,
			dimnames=list(NULL, "(Intercept)"))
}

	thelincombs=list()	
lincombMat[lincombMat==0] = NA


	noMissing=which(!someMissing)
	for(D in noMissing) {
		thisrow = lincombMat[D,]
		thisrow = as.list(thisrow[!is.na(thisrow)])
		
		thelincombs[[D]] = 
			do.call(inla.make.lincomb, thisrow)$lc

		thelincombs[[D]][[length(thelincombs[[D]])+1]] =
				list(space=list(idx=D, weight=1))
	}
	for(D in which(someMissing)) {
		thelincombs[[D]] = list(list("(Intercept)"=list(weight=1)))
	}
	
	names(thelincombs) = paste("c", 1:length(thelincombs),sep="")


	# call inla
	inlaResult = inla(formula, data=data, 
			lincomb=thelincombs, 
		#	family="poisson")
	...
		)
		
 
		

	# parameter priors for result
	
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
	params$range$prior=cbind(
			x=rangeSeq,
			y=dgamma(rateSeq, shape=ratePrior["shape"], 
					rate=ratePrior["rate"]) * (rangeSeq)^(-2) * xres(cells)
		)
	if(F) {
		plot(params$range$prior, type='l')
	
		sum(params$range$prior[,"y"])*range(diff(params$range$prior[,"x"]))
	}	
	
	precLim = 	qgamma(c(0.999,0.001), 
			shape=precPrior["shape"], 
			rate=precPrior["rate"])
	sdLim = 1/sqrt(precLim)
	sdSeq = seq(min(sdLim), max(sdLim), len=1000)
	precSeq = sdSeq^(-2)
	params$sd$prior=cbind(
			x=sdSeq,
			y=dgamma(precSeq, shape=precPrior["shape"], 
					rate=precPrior["rate"]) *2* (precSeq)^(3/2) 
	)
	if(F) {
		plot(params$sd$prior, type='l')
		
		sum(params$sd$prior[,"y"])*range(diff(params$sd$prior[,"x"]))
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
		fittedMat = as.matrix(inlaResult$summary.lincomb.derived[,-1])
		fittedMat[someMissing,] = NA
		resRasterFitted[] = fittedMat
		names(resRasterFitted) = paste("predict.", names(resRasterFitted),sep="")

	# mask	
		

	# posterior distributions
params$sd$posterior=inlaResult$marginals.hyperpar[["Precision for space"]]
params$sd$posterior[,"y"] = params$sd$posterior[,"y"] * 2*  
		params$sd$posterior[,"x"]^(3/2) 
params$sd$posterior[,"x"] = 1/sqrt(params$sd$posterior[,"x"])  
params$sd$posterior = params$sd$posterior[seq(dim(params$sd$posterior)[1],1),]		

# sum(c(0,diff(params$sd$posterior[,"x"])) * params$sd$posterior[,"y"])
# sum(c(0,diff(params$sd$prior[,"x"])) * params$sd$prior[,"y"])

params$range$posterior=inlaResult$marginals.hyperpar[["Range for space"]]
params$range$posterior[,"x"] =  xres(cells) /(params$range$posterior[,"x"])
params$range$posterior[,"y"] = params$range$posterior[,"y"] *   
		params$range$posterior[,"x"]^(-2) * xres(cells)
params$range$posterior = params$range$posterior[seq(dim(params$range$posterior)[1],1),]		

# sum(c(0,diff(params$range$posterior[,"x"])) * params$range$posterior[,"y"])
# sum(c(0,diff(params$range$prior[,"x"])) * params$range$prior[,"y"])


params$summary = inlaResult$summary.fixed
params$summary = rbind(params$summary,
		sd=c(NA, NA, 
				1/sqrt(inlaResult$summary.hyperpar["Precision for space",
								paste(c("0.975", "0.5","0.025"), "quant", sep="")
								]), 
								NA),
				range=c(NA, NA, 
						xres(cells)/inlaResult$summary.hyperpar["Range for space",
										paste(c("0.975", "0.5","0.025"), "quant", sep="")
		], 
						NA)
		)
		
params$summary["range","mean"] =xres(cells)*sum(
  1/(inlaResult$marginals.hyperpar[["Range for space"]][,"x"])*
	c(0,diff(inlaResult$marginals.hyperpar[["Range for space"]][,"x"]))*
	inlaResult$marginals.hyperpar[["Range for space"]][,"y"]
)
params$summary["sd","mean"] =sum(
		1/sqrt(inlaResult$marginals.hyperpar[["Precision for space"]][,"x"])*
				c(0,diff(inlaResult$marginals.hyperpar[["Precision for space"]][,"x"]))*
				inlaResult$marginals.hyperpar[["Precision for space"]][,"y"]
)

precGauName = "Precision for the Gaussian observations"
if(precGauName %in% names(inlaResult$marginals.hyperpar)) {
	params$summary = rbind(params$summary,
			sdNugget=c(NA, NA, 
					1/sqrt(inlaResult$summary.hyperpar[precGauName,
									paste(c("0.975", "0.5","0.025"), "quant", sep="")
							]), 
					NA)
			)
	params$summary["sdNugget","mean"] =sum(
				1/sqrt(inlaResult$marginals.hyperpar[[precGauName]][,"x"])*
						c(0,diff(inlaResult$marginals.hyperpar[[precGauName]][,"x"]))*
						inlaResult$marginals.hyperpar[[precGauName]][,"y"]
	)
	
	params$sdNugget = list(posterior=
					inlaResult$marginals.hyperpar[[precGauName]]
	)
	params$sdNugget$posterior[,"y"] = params$sdNugget$posterior[,"y"] * 2*  
			params$sdNugget$posterior[,"x"]^(3/2) 
	params$sdNugget$posterior[,"x"] = 1/sqrt(params$sdNugget$posterior[,"x"])  
	params$sdNugget$posterior = 
			params$sdNugget$posterior[seq(dim(params$sdNugget$posterior)[1],1),]		
	
	
}
		


	result=list(inla=inlaResult,
					raster=stack(resRasterRandom, resRasterFitted),
					parameters=params
			)

	result
	
}