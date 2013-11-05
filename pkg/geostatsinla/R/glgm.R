glgm = function(data,  cells, covariates=NULL, formula=NULL, 
		priorCI=NULL, shape=1, buffer=0,
		mesh=FALSE,...) {

	# list of additional arguments
	thedots = list(...)
	
	# create raster for prediction
	if(!length(grep("^Raster",class(cells)))) { 
		# cells must be an integer
		cells = as.integer(cells)
		thebbox = data@bbox
		if(buffer) {
			smallBbox = thebbox
			thebbox = thebbox + buffer*cbind(-c(1,1),c(1,1))
		} else {
			smallBbox = thebbox
		}
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
		smallBbox = bbox(cells)
	}
 
	
 	if(cells@nrows * cells@ncols > 10^6) warning("there are lots of cells in the prediction raster,\n this might take a very long time")
	

	
	# the formula
	# get rid of special character is names of data
	if(!is.null(names(data)))
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


	allterms = formulaRhs(formula,char=TRUE)
	allterms = unlist(strsplit(allterms, "\\+"))
	allterms = gsub("[[:space:]]", "", allterms)
	allterms = allterms[allterms != "1"]
	
	# get rid of offset
	theOffset = grep("^offset\\(", allterms)
	if(length(theOffset)) allterms = allterms[-theOffset]
	# get rid of other random effects
	theOtherRE = grep("^f\\(", allterms)
	if(length(theOtherRE)) allterms = allterms[-theOtherRE]
	alltermsWithF = gsub("\\)$", "", allterms)
	allterms = gsub("^factor\\(", "", alltermsWithF)
	
	
	# convert covariates to raster stack with same resulution of prediction raster.
	if(!is.null(covariates)){
		# find out which variables are factors
		theFactors = grep("^factor", alltermsWithF, value=T)
		theFactors = gsub("^factor\\(", "", theFactors)
		
		# check to see if any rasters are stored as factors
		notFactors = allterms[!allterms %in% theFactors]
		for(D in notFactors) {
			if(any(slotNames(covariates[[D]])=="data")) {
				if(any(slotNames(covariates[[D]]@data)=="isfactor")) {
					if(covariates[[D]]@data@isfactor)					
						theFactors = c(theFactors, D)
				}
			}
		}
		
		method = rep("bilinear", length(covariates))
		names(method)=names(covariates)
		
		if(!all(theFactors %in% c(names(covariates), names(data)))) {
			warning("some covariates in the model aren't in the data")
		}
		method[names(method) %in% theFactors] = "ngb" 
		
		covariatesOrig = covariates
		
		covariates = stackRasterList(covariates, cells, method=method)

		# see if any factor variables aren't coded as factors in `covariates' raster
		if(length(levels(covariates))) {
			haveLevels = as.logical(unlist(lapply(levels(covariates), length)))
			haveLevels = names(covariates)[haveLevels]
		} else {
			haveLevels = NULL
		}
		needLevels = theFactors[! theFactors%in%haveLevels]
		for(D in needLevels) {
			stuff = covariates[[D]]
			theunique = unique(stuff)
			levels(stuff) = list(data.frame(ID=theunique, CLASSNAMES=as.character(theunique)))
			covariates = stack(covariates[[-which(names(covariates)==D)]], stuff)
		}
	} else { #all covariates should be in data
		
		theFactors=NULL
	}

	
	
	# cell ID's for INLA
	cellsInla = cells
	values(cellsInla ) =  
			c(t(matrix(seq(1,ncell(cellsInla)), 
									nrow=nrow(cellsInla), ncol=ncol(cellsInla))))
	names(cellsInla) = "inlaCells"

	# create data frame for inla
	# if data is a raster
	if(length(grep("^Raster", class(data)))) {
 
	data = stack(data, covariates)
	data = stack(data, cellsInla)

	thenames = names(data)
	data=as.data.frame(data)
	names(data) = thenames
		data$space = data$inlaCells
		
		# get rid of any rows with -Inf or NA, 
		# usually offset of zero on the natural scale
		isMinusInf = apply(data, 1, function(qq) any(qq==-Inf | is.na(qq) ))
		data = data[!isMinusInf,]
		
	}  else { # data is a SpatialPointsDataFrame
		
		notInData = allterms[! allterms %in% names(data)]
		
		if(length(notInData)) {
		if(! all(notInData %in% names(covariates)))
			warning("some terms in the model are missing from both the data and the covariates")
		}	
		for(D in notInData) {
			data[[D]] = extract(covariatesOrig[[D]], data)
			# check for factors

		
		theNA = which(is.na(data[[D]]))
		
		if(any(D == theFactors)){
				data[[D]] = factor(data[[D]])
				if(length(theNA)) {
					warning("NA's in covariate ", D, " element ",
							paste(theNA, collapse=","), "\n replacing with baseline")
					data[[D]][theNA]=levels(data[[D]])[1]
				}
				
			} else {
			
				if(length(theNA)) {
					warning("NA's in covariate ", D, " elements ",
						paste(theNA, collapse=","), "\n replacing with zero")
					data[[D]][theNA]=0
			}
		}
			
				
		}	
			
 	
	 
		data$space = extract(cellsInla, data ) 
		data = data@data
	}
	# data is now a data frame.
 
	# priors for spatial standard deviation and nugget std dev.
	sdNames = unique(c("sd",grep("^sd", names(priorCI), value=TRUE)))
	# if model is Gaussian, look for prior for sdNugget
	if(!any(names(thedots)=="family")) {
		thedots$family =  "gassian"
	}
	if(thedots$family=="gaussian") {
		sdNames = unique(c(sdNames, "sdNugget"))
	}
	
	precPrior=list()
	for(Dsd in sdNames) {
			
		if(any(names(priorCI)==Dsd)) {
		obj1 = sort(priorCI[[Dsd]]^-2)
		cifun = function(pars) {
				theci = 	pgamma(obj1, shape=pars[1], 
						rate=pars[2],log.p=T)
				
				(log(0.025) - theci[1])^2 +
				(2*(log(0.975) - theci[2]))^2		

			}
	
			
			
		precPrior2=optim(c(.5,.5/mean(obj1)), cifun, 
				lower=c(0.000001,0.0000001),method="L-BFGS-B")
		names(precPrior2$par) = c("shape","rate")
		precPrior[[Dsd]] = precPrior2$par 
				
 		#pgamma(obj1, shape= precPrior["shape"], rate=precPrior["rate"],log.p=F)
		#pgamma(obj1, shape= precPrior["shape"], rate=precPrior["rate"],log.p=T)
		#log(c(0.025, 0.975))
		#precPrior2
		#pgamma(obj1, shape=precPrior["shape"], rate=precPrior["rate"],log.p=T)
		#log(c(0.025, 0.975)) 
 		#1/sqrt(qgamma(c(0.975,0.025), shape=precPrior["shape"], rate=precPrior["rate"]))
		#priorCI$sd
		
		} else {
			precPrior[[Dsd]] = c(shape=0.01, rate=0.01)
		}
	}
		
	if("range" %in% names(priorCI)) {
		if(priorCI$range[1] < xres(cells)/4) {
			priorCI$range[1] = xres(cells)/4
			warning("lower bound of range CI too small, setting it to 1/4 cell size")
			
		}
		
		# rang parameter, in terms of cells, not km.
		obj1=sort(priorCI$range/xres(cells))
		
		cifun = function(pars) {
			theci = 		pgamma(obj1, shape=pars[1], rate=pars[2], log.p=T)
			
			(theci[1] - log(0.025))^2 +
					(theci[2] - log(0.925))^2 
		}
			
		ratePrior2=optim(c(2,2/mean(obj1)), cifun, 
				lower=c(0.001,0.001),method="L-BFGS-B")
		ratePrior = ratePrior2$par
		names(ratePrior ) = c("shape","rate")
if(FALSE) {
		ratePrior
		pgamma(obj1, shape=ratePrior["shape"], rate=ratePrior["rate"])
		 qgamma(c(0.975,0.025), shape=ratePrior["shape"], rate=ratePrior["rate"])
		xres(cells)*qgamma(c(0.025,0.975), shape=ratePrior["shape"], rate=ratePrior["rate"])
	}
	} else {
		ratePrior = c(shape=0.01, rate=0.01)
	}
	spaceFormula = paste(".~.+ f(space, model='matern2d', ",
				"nrow=", nrow(cells), 
				", ncol=", ncol(cells),
				", nu=", shape, 
				", hyper = list(",
				 "range=list( param=c(",
				      paste(ratePrior, collapse=","),
				"), prior='loggamma'),",
				"prec=list( param=c(",
				paste(precPrior$sd, collapse=","),
				"),prior='loggamma')",
				" ) )" 
			)
	
	formula = update.formula(formula,	as.formula(spaceFormula))
	
	
	# create linear combinations object for prediction.
	# create formula, strip out left variable and f(...) terms
 	formulaForLincombs = unlist(strsplit(as.character(formula), "~"))
	formulaForLincombs = formulaForLincombs[length(formulaForLincombs)]
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

	

	if(!is.null(covariates)) {
		covForLincomb = stack(cellsInla, covariates)
	} else{
		covForLincomb = cellsInla
	}

	

	if(!is.null(smallBbox)) {
		covForLincomb =
				crop(covForLincomb, extent(smallBbox))
	}
	lincombCells = covForLincomb[['inlaCells']]
	names(lincombCells) = "lincombCells"
	values(lincombCells) = seq(1, ncell(lincombCells))
	covForLincomb = stack(covForLincomb, lincombCells)

	thelincombs = matrix(values(covForLincomb), 
			nrow=ncell(covForLincomb), ncol=nlayers(covForLincomb),
			dimnames = list(values(covForLincomb[['lincombCells']]),
					names(covForLincomb) ) )
	
	thelincombs = na.omit(thelincombs)
	thelincombs = as.data.frame(thelincombs)
	
	
	# if we have covariates
	if(nchar(formulaForLincombs)) {
		formulaForLincombs=as.formula(paste("~", formulaForLincombs))
	
	
	# variables in the model but not in prediction rasters
	thevars = rownames(attributes(terms(formulaForLincombs))$factors)
	varsInPredict = thevars[thevars %in% names(thelincombs)]
	cantPredict = thevars[! thevars %in% names(thelincombs)]
	theFactors2 = grep("factor\\([[:print:]]*\\)", cantPredict)
	if(length(theFactors2)) {
		temp = cantPredict
		cantPredict = cantPredict[-theFactors2]
		theFactorsInFormula = temp[theFactors2]
	}
	thelincombs[,cantPredict]= 0

	# check for factors	
 
	# convert raster data into factors using same levels as in points data
	for(D in theFactors) {
			# convert columns in lincombs to factors
			# get rid of levels not present in the observed data

			thelincombs[,D] = factor(thelincombs[,D],
					levels=levels(data[,D]))
	}
	thelincombs = na.omit(thelincombs)
	
	lincombMat = model.matrix(formulaForLincombs, thelincombs, na.action=NULL)
	
} else { # no covariates
	lincombMat = matrix(rep(1,ncell(covForLincomb)), ncol=1)
			colnames(lincombMat) = "(Intercept)"
# for some reason can't give name (Intercept) in the data.frame call.
}

lincombMatCells =  thelincombs[,c("inlaCells", "lincombCells")]

thelincombs=list()	
lincombMat[lincombMat==0] = NA

		
for(D in 1:nrow(lincombMat)) {
		thisrow = lincombMat[D, ]
		thisrow = as.list(thisrow[!is.na(thisrow)  ])
		
		thelincombs[[D]] = 
			do.call(INLA::inla.make.lincomb, thisrow)$lc

		thelincombs[[D]][[length(thelincombs[[D]])+1]] =
				list(space=list(idx=lincombMatCells[D,"inlaCells"], 
								weight=1))
}

	
	names(thelincombs) = paste("c", lincombMatCells[,"lincombCells"],sep="")

#	forInla = list(formula=formula, data=data, lincomb=thelincombs)
#	forInla = c(forInla, list(...))

#	if(any(names(forInla)=="Ntrials")) 	{
#		forInla$Ntrials = forInla$Ntrials[]
#	}
	# call inla
	
	if(!is.null(thedots$lincomb))
		thelincombs = c(thelincombs, thedots$lincomb)

	forInla = list(formula=formula, data=data, 
			lincomb=thelincombs)
	forInla = c(forInla, thedots)
	
	# if model is gaussian, add prior for nugget
	if(!is.null(precPrior$sdNugget)) {
		forInla$control.family$hyper$prec =
				list(prior="loggamma",
						param=precPrior$sdNugget
				) 
	}

 
	
	inlaResult = do.call(INLA::inla, forInla) #(formula, data=data, 
#			lincomb=thelincombs, 
 	 #	family="poisson")
	#	family="binomial",verbose=T, Ntrials=Ntrials)
#	... )


 	# parameter priors for result


	params = list(
			range = list(userPriorCI=priorCI$range, 
					priorCI = 
							xres(cells)*
							qgamma(c(0.025,0.975), 
									shape=ratePrior["shape"], 
									rate=ratePrior["rate"]),
					priorCIcells = 
							qgamma(c(0.975,0.025), 
									shape=ratePrior["shape"], 
									rate=ratePrior["rate"]),
					params.intern = ratePrior))
	rangeLim = 	qgamma(c(0.001,0.999), 
			shape=ratePrior["shape"], 
			rate=ratePrior["rate"])
	rangeLim = xres(cells) *rangeLim
	rangeSeq = seq(min(rangeLim), max(rangeLim), len=1000)
	rangeSeqCells = rangeSeq/xres(cells)
	params$range$prior=cbind(
			x=rangeSeq,
			y=dgamma(rangeSeqCells, shape=ratePrior["shape"], 
					rate=ratePrior["rate"])  / xres(cells)
	)
	if(FALSE) {
		plot(params$range$prior, type='l')
		
		sum(params$range$prior[,"y"])*range(diff(params$range$prior[,"x"]))
	}	
	

	
	for(Dsd in names(precPrior)) {
		params[[Dsd]] = list(userPriorCI=priorCI[[Dsd]], 
			priorCI = 1/sqrt(
				qgamma(c(0.975,0.025), 
						shape=precPrior[[Dsd]]["shape"], 
						rate=precPrior[[Dsd]]["rate"])),
					params.intern=precPrior[[Dsd]])
	
	precLim = 	qgamma(c(0.999,0.001), 
			shape=precPrior[[Dsd]]["shape"], 
			rate=precPrior[[Dsd]]["rate"])
	sdLim = 1/sqrt(precLim)
	sdSeq = seq(min(sdLim), max(sdLim), len=1000)
	precSeq = sdSeq^(-2)
	params[[Dsd]]$prior=cbind(
			x=sdSeq,
			y=dgamma(precSeq, shape=precPrior[[Dsd]]["shape"], 
					rate=precPrior[[Dsd]]["rate"]) *2* (precSeq)^(3/2) 
	)

	}


	
	# random into raster
	if("summary.random" %in% names(inlaResult)) {
		forRast = 	as.matrix(inlaResult$summary.random[["space"]])#[,-1])
		forRastArray = array(forRast, 
				c(nrow(cells), ncol(cells),
						dim(forRast)[2]))
		dimnames(forRastArray)[[3]] = 
				colnames(forRast)
		forRastArray = aperm(forRastArray, c(2,1,3))

		resRasterRandom = 
				brick(cells,
						nl=dim(forRastArray)[3])
		names(resRasterRandom) = paste("random.", colnames(forRast),sep="")
		
		values(resRasterRandom) = forRastArray
		
		if(buffer) {
			resRasterRandom =  
				crop(resRasterRandom, extent(smallBbox))
		# remove boundary cells from inla marginals 
		cellsSmall = crop(cellsInla, extent(smallBbox))
		cellIdSmall = values(cellsSmall)
		inlaResult$marginals.random$space = 
				inlaResult$marginals.random$space[cellIdSmall]	
		}
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



	linc = inlaResult$summary.lincomb.derived
	linc$cell = as.integer(gsub("^c", "", rownames(linc)))
	
	missingCells = which(! seq(1, ncell(lincombCells)) %in% linc$cell)
	toadd = matrix(NA, length(missingCells), dim(linc)[2], 
					dimnames=list(NULL, colnames(linc)))
	toadd[,"cell"] = missingCells
	linc = rbind(linc, toadd)
	linc = linc[order(linc[,"cell"]),]
			
#	linc = linc[,!colnames(linc)%in% c("ID","cell")]
	

resRasterFitted = 
			brick(lincombCells,
					nl=dim(linc)[2])

	values(resRasterFitted) = as.matrix(linc)
	names(resRasterFitted) = paste("predict.", colnames(linc),sep="")

	# Add in empty lists for the marginals of missing cells
	missingNames = paste("c", missingCells, sep="")
	missingMarginals = vector("list", length(missingNames))
	names(missingMarginals) = missingNames

	
	inlaResult$marginals.lincomb.derived = c(	
			inlaResult$marginals.lincomb.derived,
			missingMarginals)
	
	inlaResult$marginals.lincomb.derived = 
			inlaResult$marginals.lincomb.derived[
					paste("c", seq(1,ncell(resRasterFitted)),sep="")
					]
			
	
	
	# posterior distributions




params$range$posterior=inlaResult$marginals.hyperpar[["Range for space"]]
params$range$posterior[,"x"] =  xres(cells) * params$range$posterior[,"x"]
params$range$posterior[,"y"] = params$range$posterior[,"y"] / xres(cells)


# sum(c(0,diff(params$range$posterior[,"x"])) * params$range$posterior[,"y"])
# sum(c(0,diff(params$range$prior[,"x"])) * params$range$prior[,"y"])


params$summary = inlaResult$summary.fixed

thecols = paste(c("0.975", "0.5","0.025"), "quant", sep="")

thesd = c(
		sdNugget= grep("^Precision[[:print:]]*Gaussian observations$", 
				names(inlaResult$marginals.hyperpar), value=TRUE),
		sd = grep("^Precision[[:print:]]*space$", 
				names(inlaResult$marginals.hyperpar), value=TRUE)
)

params$summary = rbind(params$summary,
		matrix(NA, nrow=length(thesd)+1, ncol=ncol(params$summary),
				dimnames = list(c("range", names(thesd)), NULL))
)


# convert precisions to standard deviations
for(Dsd in names(thesd)) {
	
	params[[Dsd]]$posterior=
			inlaResult$marginals.hyperpar[[thesd[Dsd]]]
	params[[Dsd]]$posterior[,"y"] = params[[Dsd]]$posterior[,"y"] * 2*  
			params[[Dsd]]$posterior[,"x"]^(3/2) 
	params[[Dsd]]$posterior[,"x"] = 1/sqrt(params[[Dsd]]$posterior[,"x"])  
	params[[Dsd]]$posterior = params[[Dsd]]$posterior[
			seq(dim(params[[Dsd]]$posterior)[1],1),]		

	params$summary[Dsd, thecols] = 
				1/sqrt(inlaResult$summary.hyperpar[
								thesd[Dsd],rev(thecols)])

		
	params$summary[Dsd,"mean"] =sum(
		1/sqrt(inlaResult$marginals.hyperpar[[thesd[Dsd]]][,"x"])*
			c(0,diff(inlaResult$marginals.hyperpar[[thesd[Dsd]]][,"x"]))*
				inlaResult$marginals.hyperpar[[thesd[Dsd]]][,"y"]
	)
}

# put range in summary, in units of distance, not numbers of cells
thecolsFull =c("mean","sd",thecols,"mode") 
params$summary["range",thecolsFull]=				
		xres(cells)*
		inlaResult$summary.hyperpar[
				"Range for space",
				thecolsFull
		]


resRaster=stack(resRasterRandom, resRasterFitted)

#	if(!is.null(smallBbox))
#		resRaster = crop(resRaster, extent(smallBbox))
		
	result=list(inla=inlaResult,
					raster=resRaster,
					parameters=params
			)

	result
	
}