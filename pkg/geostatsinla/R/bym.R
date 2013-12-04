`inla.nb.to.graph` = function(adjMat, graph.file="graph.dat")
{
	## A function for converting GeoBUGS adjacency data into the INLA
	## graph format. Kindly provided by Aki Havunlinna tkk.fi; thanks.
	
	fd = file(graph.file,  "w")
	len <- length(adjMat)

	cat(len, '\n', file=fd)
	k = 1L
	for(i in 1L:len) {
		num=adjMat[[i]]
		lnum = length(num)
		if (num[1] == "0" | lnum==0 ) {
			cat(i, "0", "\n", file=fd)
		} else {
			cat(i, lnum, 
					num, "\n", file = fd)
		}
	}
	close(fd)

	
	region.index = 1:len
	region.id = attributes(adjMat)$region.id
	if(is.null(region.id))
		region.id = region.index
	names(region.index) = as.character(region.id)
		
	return(region.index)
}

bym = function(data, ...) {
	UseMethod("bym")
	
}


bym.SpatialPolygonsDataFrame = function(data, 
		formula=observed ~ offset(logExpected), 
		priorCI=list(sdSpatial = c(0.01, 10), sdIndep = c(0.01, 10)), 
		family="poisson",
		region.id,adjMat, formula.fitted=formula,
		...
	) {
	
	if(missing(region.id)) {
		region.id="region.id"
		data[[region.id]] = seq(1,dim(data)[1])
	}
		
	if(missing(adjMat))
		adjMat=spdep::poly2nb(data, row.names =  data[[region.id]] )
 
	result = bym.data.frame(data=data.frame(data),
			formula=formula, priorCI=priorCI, family=family,
			region.id=region.id, adjMat=adjMat,formula.fitted=formula.fitted
		, ...
			)
	if(any(names(result)=="data")) {		
 	# merge data back into SPDF
		data@data = result$data[as.character(data[[region.id]]),]
		result$data = data
	}
	
	result
}

bym.data.frame = function(data,   formula, 
		priorCI=list(sdSpatial = c(0.01, 10), sdIndep = c(0.01, 10)), 
		family="poisson",
		region.id, adjMat,
		formula.fitted=formula,
		...
	) {
		#neighbourhood structure
		if(missing(region.id)) {
			region.id="region.id"
			data[[region.id]] = seq(1,dim(data)[1])
		}
		
		graphfile=tempfile()
		
		# if using windows, replace back slashes with forward slashes...
		graphfile = gsub("\\\\", "/", graphfile)
		
		region.index = inla.nb.to.graph(adjMat, graphfile)

		# check for data regions missing from adj mat
		data[[region.id]] = as.character(data[[region.id]])
		if(!all(data[[region.id]] %in% names(region.index)))
			warning("regions in data missing from adjacency matrix")
		data$region.indexS = data$region.indexI = region.index[data[[region.id]]]
		
		# priors
		if(!all(c("sdSpatial","sdIndep")%in% names(priorCI))) {
			warning("priorCI needs elements sdSpatial and sdIndep")
		}
		
		cifun = function(pars) {
			theci = 	pgamma(obj1, shape=pars[1], rate=pars[2],log.p=T)
			
			(log(0.025) - theci[1])^2 +
					(2*(log(0.975) - theci[2]))^2		
			
		}

		precPrior = list()
		for(D in c("sdSpatial","sdIndep")) {
			obj1 = sort(priorCI[[D]]^-2)
		
			startMean = mean(obj1)
			startSD = diff(obj1)/4
			startShape = startSD^2/startMean^2
			
			precPrior2=optim(c(startShape,startShape/startMean), cifun, 
					lower=c(0.000001,0.0000001),method="L-BFGS-B",
					control=list(parscale=c(startShape,startShape/startMean)))
			precPrior[[D]] = precPrior2$par
			names(precPrior[[D]] ) = c("shape","rate")
			
			pgamma(obj1, shape= precPrior[[D]]["shape"], rate=precPrior[[D]]["rate"],log.p=F)
			pgamma(obj1, shape= precPrior[[D]]["shape"], rate=precPrior[[D]]["rate"],log.p=T)
			log(c(0.025, 0.975))
			precPrior2
			pgamma(obj1, shape=precPrior[[D]]["shape"], rate=precPrior[[D]]["rate"],log.p=T)
			log(c(0.025, 0.975)) 
			1/sqrt(qgamma(c(0.975,0.025), shape=precPrior[[D]]["shape"], rate=precPrior[[D]]["rate"]))
			priorCI[[D]]
			
		} 
		
		#f(CSDUID, model = "bym", graph.file = "nb.graph", 
		#		+     param = c(prior.iid, prior.besag), values = CSDUID)
		
		
		bymTerm = paste(
				".~.+f(region.indexS, model='besag', graph='",
				graphfile,
				"', hyper = list(theta=list(param=c(",
				paste(precPrior[["sdSpatial"]], collapse=","), ")))) ",
				"+ f(region.indexI, model='iid',",
				"hyper=list(theta=list(param=c(", 		
				paste(precPrior[["sdIndep"]], collapse=","),
				"))))",
			sep="")
		
		formula = update(formula, as.formula(bymTerm))

		# linear combinations


	#check to see if some regions don't have data.  If so they'll have to be added so 
	# the independent random effect will be computed
	notInData = region.index[!region.index %in% data$region.indexI]
	if(length(notInData)) {
		dataToAdd = data[rep(1,length(notInData)),]
		rownames(dataToAdd) = paste("missing",notInData,sep="")
		dataToAdd[,"region.indexI"] = dataToAdd[,"region.indexS"]=notInData
	# set response to missing
		dataToAdd[,formulaLhs(formula)] = NA		
		data = rbind(data, dataToAdd)
	}

	#	theDup = duplicated(data$region.indexS)
		
		inlaLincombs = list()
		# random effects
		for(D in 1:length(region.index)) {
			inlaLincombs[[D]] = list(
					list(region.indexS=
									list(idx=region.index[D], weight=1)),
					list(region.indexI = 
									list(idx=region.index[D], weight=1))
							)
		}
		names(inlaLincombs) = paste("bym", names(region.index),sep="_")

		
		# fitted values
 
formulaForLincombs = formulaRhs(formula.fitted,char=TRUE)
 
# get rid of f(stuff) in formula
formulaForLincombs =
		gsub("f\\([[:print:]]*\\)", "", formulaForLincombs)
# get rid of offset(stuff)
formulaForLincombs =
		gsub("offset\\([[:print:]]*\\)[[:space:]]?($|\\+)", "", formulaForLincombs)

# convert multiple + to a single +
formulaForLincombs = gsub(
		"\\+[[:space:]]?\\+([[:space:]]?\\+)?", "+",
		formulaForLincombs)
# strip out trailing or leading +
formulaForLincombs = gsub("\\+[[:space:]]?$|^[[:space:]]?\\+[[:space:]]+", "", formulaForLincombs)


startIndex = length(region.index)
 
	if(nchar(formulaForLincombs) & formulaForLincombs != "1" &
			!length(grep("^[[:space:]]+$", formulaForLincombs))) {
 
		formulaForLincombs=as.formula(
			paste("~", paste(c("1",formulaForLincombs),collapse="+"))
		)
 
		# remove regions in the data set twice
		theDuplicated = duplicated(data$region.indexS)
		notDuplicated = which(!theDuplicated)
		
		# reorder the matrix by region ID
		dataOrder = data[notDuplicated,]
		dataOrder = dataOrder[!dataOrder$region.indexI %in% notInData,]
		dataOrder = dataOrder[order(dataOrder$region.indexI),]

		
		lincombFrame = model.frame(formulaForLincombs, dataOrder,
				na.action=na.omit)

		SregionFitted = dataOrder[rownames(lincombFrame),"region.indexI"]
		names(SregionFitted) = dataOrder[rownames(lincombFrame),region.id]
		
		
		lincombMat = model.matrix(formulaForLincombs, lincombFrame)
		
		lincombMat[lincombMat==0]= NA
		thelincombs = inla.make.lincombs(as.data.frame(lincombMat))
		for(D in seq(1,length(SregionFitted))) {	
			inlaLincombs[[D+startIndex]] = 
				c(
					list(list(region.indexS=
								list(idx=SregionFitted[D], weight=1))),
					list(list(region.indexI= 
								list(idx=SregionFitted[D], weight=1))) 
				)
			if(length(thelincombs)>=D)
				inlaLincombs[[D+startIndex]] = c(thelincombs[[D]],
					inlaLincombs[[D+startIndex]])
		}
		names(inlaLincombs)[seq(startIndex+1, len=length(SregionFitted))] =
				paste("fitted",names(SregionFitted),sep="_")
	} else { # add only intercept to predictions
		formulaForLincombs = ~1
		lincombMat = data.frame(x=rep(1,length(region.index)))
		SregionFitted = region.index
		for(D in 1:length(region.index)) {	
			inlaLincombs[[D+startIndex]] = 
					 list(
							list("(Intercept)" = list(weight=1)),
							list(region.indexS=
											list(idx=region.index[D], weight=1)),
							list(region.indexI = 
											list(idx=region.index[D], weight=1))
					)
		}
		names(inlaLincombs)[seq(startIndex+1, len=length(region.index))] =
				paste("fitted",names(region.index),sep="_")
	}
 



	# run inla!		
	inlaRes = inla(formula, data=data , family=family,
			lincomb=inlaLincombs, ...)
 	
	if(all(names(inlaRes)=="logfile"))
		return(c(list(formula=formula, data=data,
						family=family, 
						lincomb=inlaLincombs, 
						ldots = list(...)),
						inlaRes)
	)
	

	# posterior distributions of random effect (spatial + independent)
	thebym = inlaRes$summary.lincomb.derived[
			grep("^bym_", rownames(inlaRes$summary.lincomb.derived)),]
 	
	inlaRes$marginals.bym = inlaRes$marginals.lincomb.derived[
			grep("^bym_", names(inlaRes$marginals.lincomb.derived), value=TRUE)
			]
 	

	thebym = thebym[,!names(thebym) %in% c("ID","kld")]
	colnames(thebym) = paste("random.",colnames(thebym),sep="")
	rownames(thebym) = gsub("^bym_", "", rownames(thebym))	
	names(inlaRes$marginals.bym) = gsub("^bym_", "", 
			names(inlaRes$marginals.bym) )
	# make sure they're in the correct order
	thebym = thebym[names(region.index),]
	inlaRes$marginals.bym = inlaRes$marginals.bym[names(region.index)] 

	
	# fitted values, some regions dont have them if covariates are missing
	theFitted = inlaRes$summary.lincomb.derived[
			grep("^fitted_", rownames(inlaRes$summary.lincomb.derived)),]
	inlaRes$marginals.fitted.bym = inlaRes$marginals.lincomb.derived[
			grep("^fitted_", names(inlaRes$marginals.lincomb.derived))
			]
			
	names(inlaRes$marginals.fitted.bym)= 
			gsub("^fitted_", "", names(inlaRes$marginals.fitted.bym))
	rownames(theFitted) = gsub("^fitted_", "", rownames(theFitted))
	
	meanExp = unlist(
					lapply(inlaRes$marginals.fitted.bym, 
							function(qq) {
								sum(
										exp(qq[,"x"])*c(0,diff(qq[,"x"]))*qq[,"y"]	
								)
							})
				) # end unlist
	meanExp[meanExp==Inf]=NA
	theFitted = cbind(theFitted, exp = meanExp[rownames(theFitted)])

	
	# E inv logit(lincombs)

	if(length(grep("binomial",inlaRes$.args$family))) {
		invlogit=unlist(
				lapply(inlaRes$marginals.fitted.bym, 
						function(qq) {
							eqqx = exp(qq[,"x"])
							sum(
									eqqx/(1+eqqx)*c(0,diff(qq[,"x"]))*qq[,"y"]	
							)
						})
		)
		theFitted = cbind(theFitted, invlogit = invlogit[rownames(theFitted)])
	}

	
	theFitted = theFitted[,!names(theFitted) %in% c("ID","kld")]
	colnames(theFitted) = paste("fitted.",colnames(theFitted),sep="")
	
	# merge fitted falue summary into BYM
	thebym = cbind(thebym, matrix(NA, dim(thebym)[1], dim(theFitted)[2],
					dimnames=list(NULL, names(theFitted))))
	thebym[rownames(theFitted), colnames(theFitted)] = theFitted
	
	
	# make fitted marginals list have same names and order as fitted radom
	notInFitted = names(inlaRes$marginals.bym) [
		!names(inlaRes$marginals.bym)%in%names(inlaRes$marginals.fitted.bym)]

	toAdd = replicate(length(notInFitted),NULL,simplify=FALSE)
	names(toAdd) = notInFitted
	inlaRes$marginals.fitted.bym = c(inlaRes$marginals.fitted.bym, 
			 toAdd)
	 inlaRes$marginals.fitted.bym = inlaRes$marginals.fitted.bym[
			 names(inlaRes$marginals.bym) 
			 ]
	
	# the parameters
	params=list()
	params$summary = inlaRes$summary.fixed
	
	quantNames = grep("quant$", colnames(params$summary), value=TRUE)
	revQuant = rev(quantNames)	
	
	for(D in c("S","I")) {
		
		Dname = grep( paste("^sd",D,sep=""),names(priorCI),value=TRUE)
		
		params[[Dname]] = list(
				userPriorCI=priorCI[[Dname]], 
				priorCI = 
						1/sqrt(
								qgamma(c(0.975,0.025), 
										shape=precPrior[[Dname]]["shape"], 
										rate=precPrior[[Dname]]["rate"])),
				params.intern=precPrior[[Dname]])
	
		
		imname = paste("Precision for region.index",D,sep="")
		params[[Dname]]$posterior=
				inlaRes$marginals.hyperpar[[
					imname	
				]]
		params[[Dname]]$posterior[,"y"] = params[[Dname]]$posterior[,"y"] * 2*  
				params[[Dname]]$posterior[,"x"]^(3/2) 
		params[[Dname]]$posterior[,"x"] = 1/sqrt(params[[Dname]]$posterior[,"x"])  
		params[[Dname]]$posterior = 
				params[[Dname]]$posterior[seq(dim(params[[Dname]]$posterior)[1],1),]		
		
		
		precLim = range( inlaRes$marginals.hyperpar[[
								paste("Precision for region.index",D,sep="")
						]][,1] ) 
		precLim = precLim * c(0.8, 1.2)		
		sdLim = 1/sqrt(precLim)
		sdSeq = seq(min(sdLim), max(sdLim), len=100)
		precSeq = sdSeq^(-2)
		params[[Dname]]$prior=cbind(
				x=sdSeq,
				y=dgamma(precSeq, shape=precPrior[[Dname]]["shape"], 
						rate=precPrior[[Dname]]["rate"]) *2* (precSeq)^(3/2) 
		)
		
		
		thesummary = inlaRes$summary.hyperpar[imname, ,drop=FALSE]
		thesummary[,quantNames] = 1/sqrt(thesummary[,revQuant])
		
		thesummary[,"mean"] =sum(
				1/sqrt(inlaRes$marginals.hyperpar[[imname]][,"x"])*
						c(0,diff(inlaRes$marginals.hyperpar[[imname]][,"x"]))*
						inlaRes$marginals.hyperpar[[imname]][,"y"]
		)
		thesummary[,"sd"] = NA
		rownames(thesummary) = Dname
 
		
		donthave = colnames(params$summary)[
				!colnames(params$summary) %in% colnames(thesummary)]

		thesummary = cbind(thesummary, matrix(NA, nrow(thesummary), length(donthave),
							dimnames=list(NULL, donthave)))
 
			
		params$summary = rbind(params$summary, 
				thesummary[,colnames(params$summary),drop=FALSE]
				)		
}
# sum(c(0,diff(params$sd$posterior[,"x"])) * params$sd$posterior[,"y"])
# sum(c(0,diff(params$sd$prior[,"x"])) * params$sd$prior[,"y"])
if(FALSE) {
	plot(params$sdSpatial$prior, type='l',lwd=4)
	lines(params$sdSpatial$posterior,col='orange',lwd=2)
	
	sum(params$sdSpatial$prior[,"y"])*range(diff(params$sdSpatial$prior[,"x"]))
}	

 



	return(list( inla=inlaRes, data=thebym, parameters=params))
}