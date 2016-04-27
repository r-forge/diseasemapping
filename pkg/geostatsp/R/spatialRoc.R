spatialRocPolyTemplate = function(
		truth, fit
) {
	
	breaksInterior = breaks[seq(2, length(breaks)-1)]

	toRasterize = fit[[1]]$data
	toRasterize$fitID = 1:length(toRasterize)
	template = rasterize(
			toRasterize,
  		truth,
			field='fitID'
			)

	names(template) = 'fitID'
	
	template
	
}

spatialRocRasterTemplate = function(
		truth, fit
) {
	
	template = raster(fit[[1]]$raster)
	values(template) = seq(1, ncell(template))
	names(template) = 'fitID'
	
	# remove cells with no predictions
	template = mask(template, fit[[1]]$raster$predict.mean)
	
	# and cells with no truth
	template = mask(template, truth[[1]])
	
	templateID = stackRasterList(
			list(fitID = template), truth, method='ngb'
	)

	templateID
}

spatialRocSims = function(
		truthCut, marginals, templateID, 
		breaks, prob=seq(0,1,len=100)
		) {
			

			breaksInterior = breaks[seq(2, length(breaks)-1)]
			
			

			Slevels = seq(1,length(breaks)-1)

			truthCdf = zonal(
					truthCut,
					templateID,
					function(x,...) {
						ecdf(x)(Slevels)
						}
					)
			rownames(truthCdf) = truthCdf[,'zone']
			truthCdf = truthCdf[,grep('zone', colnames(truthCdf), invert=TRUE)]
			colnames(truthCdf) = paste(
					rep(names(truthCut), 
							rep(length(Slevels),nlayers(truthCut))
							), 
					rep(Slevels, nlayers(truthCut)),
					sep="_"
			)
			truthCdf = array(
					truthCdf,
					c(nrow(truthCdf), length(Slevels), nlayers(truthCut)),
					dimnames = list(
							rownames(truthCdf),
							Slevels,
							names(truthCut)
							)
					)
			truthCdfUpper = 1-truthCdf
			
			idTable = table(values(templateID))		
			idMatrix = matrix(
					idTable[dimnames(truthCdf)[[1]]],
					length(idTable),
					length(prob)
					)
			
			truthCells = array(
					idTable[dimnames(truthCdf)[[1]]],
					dim=dim(truthCdf)
			)
			
			truthOver = truthCdfUpper * truthCells 
			truthUnder = truthCdf * truthCells 
			

				
			tP = tN = array(NA,
					c(length(prob), length(breaksInterior), length(fit)),
					dimnames= list(
							prob,breaksInterior, names(fit)
							)
					)		
					
					
			for(Dsim in 1:length(fit)) {
				

				notNull = which(!unlist(lapply(
										marginals[[Dsim]], is.null
										)))
				
				pMat = simplify2array(
						lapply(
						marginals[[Dsim]][notNull], 
						INLA::inla.pmarginal, 
						q=breaksInterior)
				)
				rownames(pMat) = as.character(breaksInterior)
				# upper tail probabilities
				pMat = pMat[,match(names(marginals[[Dsim]]), colnames(pMat))]
				pMat = t(pMat)
				
				
				pMat = pMat[as.numeric(dimnames(truthCdf)[[1]]),]
				pMat = 1-pMat

				for(Dlevel in 1:ncol(pMat)) {

					predOver = outer(pMat[,Dlevel], prob, '>')
				
					tP[,Dlevel,Dsim] = apply(predOver*truthOver[,Dlevel,Dsim],2,sum, na.rm=TRUE)
					tN[,Dlevel,Dsim] = apply((1-predOver)*truthUnder[,Dlevel,Dsim],2,sum, na.rm=TRUE)
					
				}
				
				
			}
				
			

			return(list(
							tP = tP, 
							allP = apply(truthOver, c(2,3), sum), 
							tN = tN, 
							allN = apply(truthUnder, c(2,3), sum)
							))
			
			
			

			
		}


spatialRoc = function(fit, 
		rr=c(1,1.2, 1.5,2), truth, 
		border=NULL, random=FALSE,
		prob = NULL, spec = seq(0,1,by=0.01)){
	
	if(is.null(prob)){
		prob = 2^seq(-1, -12, len=25)
		prob = sort(unique(c(prob, 1-prob)))
	}
	
	if(any(names(fit)=='inla')){
		fit = list(fit)
	}
	
	
	breaks = c(-Inf, log(rr), Inf)
	
	if('raster'%in% names(truth))
		truth = truth$raster
	
	if(random) {
		truthVariable = 'random'
		truth = truth[[
				intersect(
						paste(truthVariable, c('',1:length(fit)), sep=''),
						names(truth)
				)
		]]
	} else { 
		truthVariable = 'relativeIntensity'
		
		truth = truth[[
				intersect(
						paste(truthVariable, c('',1:length(fit)), sep=''),
						names(truth)
				)
		]]
		tname = names(truth)
		truth = log(truth)
		names(truth) = tname
	} 

	if(is.null(border))
		border = fit[[1]]$data
		
	if(!is.null(border))
		truth = mask(truth, border)
	

	truthCut = cut(truth, breaks=breaks)
	names(truthCut) = names(truth)
	
	
	
	if(any(names(fit[[1]]) =='raster')){
		# lgcp or glgm
		
		templateID = spatialRocRasterTemplate(
			truthCut, fit
		) 
	
	
		if(random) {
			marginals = lapply(
				fit, function(x){
					 x$inla$marginals.random$space
				}
				)
		} else {
			marginals = lapply(
				fit, function(x){
					x$inla$marginals.predict
				}
			)
		}
		
	} else { # bym model

		templateID = spatialRocPolyTemplate(
				truthCut, fit
		) 
		
		if(random) {
			marginals = lapply(
					fit, function(x){
					  x$inla$marginals.bym
					}
			)
		} else {
			marginals = lapply(
					fit, function(x){
					  x$inla$marginals.fitted.bym
					}
			)
		}
		
	}

	res = spatialRocSims(
			truthCut, marginals, templateID, 
			breaks, prob
		)	
	
	
	allParray = array(
					res$allP,
					dim = c(dim(res$allP), dim(res$tP)[1]),
					dimnames = c(
							dimnames(res$allP),
							dimnames(res$tP)[1]
							)
			)

	allNarray = array(
			res$allN,
			dim = dim(allParray),
			dimnames = dimnames(allParray)
	)

	allParray = aperm(allParray, c(3,1,2))		
	allNarray = aperm(allNarray, c(3,1,2))		
	
	

			bySim = abind::abind(
			spec = res$tP / allParray[,-dim(allParray)[2], ,drop=FALSE],
			sens = res$tN / allNarray[,-dim(allNarray)[2],,drop=FALSE],
			along=4
			)

			result = apply(bySim, c(1,2,4), mean)

			dimnames(result)[[2]] = rr
					
			
			if(!is.null(spec)) {
				
				resultOut = matrix(
						NA,
						length(spec), 
						dim(result)[2]+1,
						dimnames = list(
								1:length(spec),
								c('spec', dimnames(result)[[2]])
								)
						)
				resultOut[,'spec'] = spec
				for(D in dimnames(result)[[2]]) {
					resultOut[,D] = approx(
							result[,D,'spec'],
							result[,D,'sens'],
							xout = spec 
							)$y
				}
				resultOrig = result
				result = resultOut
				attributes(result)$orig = resultOrig
			}
			
			attributes(result)$sim = bySim
	
  result
}
