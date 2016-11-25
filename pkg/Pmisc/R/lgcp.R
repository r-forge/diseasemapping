#' Convert lgcp results to geostatsp format
#' 
#' @description Results from \code{\link[lgcp]{lgcp}} are converted 
#' to a format comparable to that produced by \code{\link[geostatsp]{lgcp}}
#' from the \code{geostatsp} package
#' 
#' @param x result from \code{\link[lgcp]{lgcp}} 
#' @param dirname folder with detailed results
#' 
#' @import raster
#' @import sp
#' 
#' @export
lgcpToGeostatsp = function(x, dirname = x$gridfunction$dirname) {
# convert lgcp results to geostatsp-type stuff
	
	bRes = x
	if(is.null(dirname)) dirname = tempdir()
	
	randomFile = file.path(bRes$gridfunction$dirname, "simout.nc")
	randomFileOut = file.path(dirname, "random.grd")
	fixedFile = file.path(dirname, 'fixed.grd')
	intensityFile = file.path(dirname, 'intensity.grd')
	summaryFile = file.path(dirname, 'summary.grd')
	maskFile = file.path(dirname, 'mask.grd')
	offsetFile = file.path(dirname, 'offset.grd')
	
	gsize = diff(bRes$ncens[1:2])/2
	scaleIntensity = sum(x$poisson.offset, na.rm=TRUE)
	if(!length(scaleIntensity))
		scaleIntensity = 1
	
	random = brick(randomFile)
	random = setMinMax(random)
	extent(random) = extent(
			c(bRes$xyt$window$xrange,
					bRes$xyt$window$yrange)
	)
	random = writeRaster(random, 
			filename = randomFileOut,
  		overwrite = file.exists(randomFileOut)
	)
	
	rasterMask = raster(
			t(bRes$cellInside[,ncol(bRes$cellInside):1]),
			template=random
	)
	rasterMask = subs(rasterMask, data.frame(1,1))
	names(rasterMask) = 'mask'
	rasterMask = writeRaster(rasterMask, maskFile,
			overwrite = file.exists(maskFile)
	)
	
	scaleOffset = sum(log(res(random)))
	

			
	smallZ = array(
			bRes$Z,
			c(bRes$M, bRes$ext, bRes$N, bRes$ext, ncol(bRes$Z)),
			dimnames = list(
					bRes$mcens,
					paste("a", 1:bRes$ext, sep="_"),
					bRes$ncens,
					paste("b", 1:bRes$ext, sep="_"), 
					names(bRes$glmfit$coefficients)
	)
	)
	smallZ = matrix(smallZ[,1,,1,], 
			ncol=ncol(bRes$Z), 
			dimnames = list(NULL, dimnames(smallZ)[[5]])
	)
	smallZ[as.vector(bRes$cellInside)==0, ] = NA
	
	fixed = tcrossprod(smallZ, bRes$betarec) + scaleOffset
	fixed = array(fixed, 
			c(bRes$M, bRes$N, nrow(bRes$betarec)))
	fixed = fixed[,dim(fixed)[2]:1,]
	fixed = aperm(fixed, c(2,1,3))
#	fixed = fixed + log(scaleIntensity)
	fixed = brick(fixed)
	extent(fixed) = extent(random)
	fixed = setMinMax(fixed)
	
	fixed = writeRaster(fixed, 
			filename = fixedFile,
  		overwrite = file.exists(fixedFile)
	)
	
	relativeIntensity = exp(fixed + random)
	relativeIntensity = setMinMax(relativeIntensity)
	
	relativeIntensity = writeRaster(
			relativeIntensity,
			filename = intensityFile,
			overwrite = file.exists(intensityFile)
	)
	
	summaryFunction = function(x, prefix = '',...) {
		Sprob=c(0.025, 0.5, 0.975)
		res = quantile(x, na.rm=TRUE, prob = Sprob)
		names(res) = paste('xx.',Sprob, "quant", sep='')
		c('xx.mean' = mean(x, na.rm=TRUE), res)
	}				
	
	randomSummary = calc(random, summaryFunction)
	names(randomSummary) = gsub("^xx", 'random', names(randomSummary))
	intensitySummary = calc(relativeIntensity, summaryFunction)
	names(intensitySummary) = gsub("^xx", 'relativeIntensity', 
			names(intensitySummary))
	
	summaryStack = stack(
			randomSummary, intensitySummary
	)
	
	offset = raster(random)
	
	values(offset) = t(x$poisson.offset[
					1:x$M, x$N:1
			])
	offset = offset*exp(-scaleOffset)
	offset = brick(
			offset, log(offset)
	)
	names(offset) = c('offsetExp','offsetLog')			
	
	summaryStack = stack(
			summaryStack, offset
	)
	summaryStack = setMinMax(summaryStack)
	
	summaryStack = mask(summaryStack, rasterMask)
	if(!all(is.na(values(summaryStack[[1]])))) {
		summaryStack = raster::trim(summaryStack,
				values = c(NA, NaN),
				filename=summaryFile,
				overwrite = file.exists(summaryFile)
		)
	}
	
	
	
	if(is.null(colnames(bRes$betarec))) {
		colnames(bRes$betarec) = names(bRes$glmfit$coef)
	}
	
	samples = list(
			relativeIntensity = relativeIntensity,
			random = random,
			fixed = fixed,
			mask=rasterMask,
			beta = bRes$betarec,
			sd = exp(bRes$etarec[,1]),
			range = 2*exp(bRes$etarec[,2])
	)
	samples$beta[,1] =
			samples$beta[,1] - samples$sd^2/2 + scaleOffset
	
	summaryTable = lapply(
			c(as.list(as.data.frame(samples$beta)), 
					samples['range'], 
					samples['sd']), 
			summaryFunction)
	
	summaryTable = t(simplify2array(summaryTable))
	colnames(summaryTable) = gsub(
			'^xx\\.', '', colnames(summaryTable)			
	)
	summaryTable = rbind(summaryTable,
			shape=c(environment(bRes$covFct)$additionalparameters, 
					rep(NA, ncol(summaryTable)-1))
	)
	
	xSeq = seq(0,
			exp(bRes$priors$etaprior$mean[1] +
							4*sqrt(bRes$priors$etaprior$variance[1,1])),
			len=1000)
	
	sd = list(
			posterior = do.call(cbind, 
					density(samples$sd)[c('x','y')]),
			prior = cbind(
					x = xSeq,
					y = stats::dlnorm(xSeq, 
							bRes$priors$etaprior$mean[1], 
							sqrt(bRes$priors$etaprior$variance[1,1])
					)
			)
	)
	
	
	xSeq = seq(0,
			exp(bRes$priors$etaprior$mean[2] +
							4*sqrt(bRes$priors$etaprior$variance[2,2])),
			len=1000)
	
	range = list(
			posterior = do.call(cbind, 
					density(samples$range)[c('x','y')]),
			prior = cbind(
					x = xSeq*2,
					y = stats::dlnorm(xSeq, 
							bRes$priors$etaprior$mean[2], 
							sqrt(bRes$priors$etaprior$variance[2,2])
					)/2
			)
	)
	
	list(
			summary = summaryTable,
			parameters = list(
					sd = sd,
					range = range),
			raster = summaryStack,
			samples = samples,
			filenames = c(
					random=randomFile,
					fixed=fixedFile,
					intensity=intensityFile,
					summary=summaryFile)
	)
}		



getZmatPopulation = function(
		events, 
		population, covariates, 
		border=population, fact=4,
		polyolay = list(
				Ncell=32
		),
		verbose=TRUE
) {
	
	if(!missing(border) & !missing(events))
		events = events[!is.na(sp::over(events, methods::as(border,'SpatialPolygons')))]
	
	if(!missing(events)) {
		eventsPPP = spatstat::as.ppp(events)
	} else {
		eventsPPP = NULL
	}
	
	if(all(class(polyolay)=='list')) {
		if(missing(border)){
			warning("border must be supplied if polyolay is a list")
		}
		if(!length(polyolay$ext))
			polyolay$ext = 2
		if(length(polyolay$Ncell)) {
			polyolay$cellwidth = min(abs(apply(bbox(border),1,diff))/polyolay$Ncell)
		}
		if(!length(polyolay$cellwidth))
			warning("polyolay must contain either cellwidth or Ncell")
		
		border2 = border
		border2$junk = 1
		border2 = border2[,'junk']
		border2@data = lgcp::assigninterp(df = border2@data, vars = "junk", 
				value = "ArealWeightedSum")
		
		if(verbose) cat("\nrunning getpolyol")
		polyolay <- lgcp::getpolyol(
				data=eventsPPP,
        regionalcovariates=border2,
        cellwidth=polyolay$cellwidth,
        ext=polyolay$ext)
		if(verbose) cat("\ndone\n")
		
	}						
	
	rasterMaskExt = raster(
			t(polyolay$gridobj$cellInside[,ncol(polyolay$gridobj$cellInside):1]),
			polyolay$gridobj$mcens[1],
			max(polyolay$gridobj$mcens),
			polyolay$gridobj$ncens[1],
			max(polyolay$gridobj$ncens)
	)
	if(!missing(population)) {
		projection(rasterMaskExt) = projection(population)			
	} else if(!missing(border)) {
		projection(rasterMaskExt) = projection(border)				
	}
	
	extentNotExt = extent(
			polyolay$gridobj$mcens[1]-1,
			mean(polyolay$gridobj$mcens[c(0,1)+polyolay$gridobj$M/polyolay$ext]),
			polyolay$gridobj$ncens[1]-1,
			mean(polyolay$gridobj$ncens[c(0,1)+polyolay$gridobj$N/polyolay$ext])
	)		
	
	rasterMask = crop(rasterMaskExt, extentNotExt)
	
	rasterFine = disaggregate(rasterMask, fact)
	
	intercept = raster(rasterMask)
	values(intercept) = 1
	names(intercept) = 'intercept'
	
	if(!missing(population)) {
		if(verbose) cat("\nrasterizing population")
		popFine = geostatsp::spdfToBrick(
				list(population),
				rasterFine,
				pattern='^expected$'
		)
		if(verbose) cat("\ndone\n")
		if(!missing(border)) {
			popFine = mask(popFine, border)		
		}
		popCoarse = aggregate(popFine, fact, fun=sum, na.rm=TRUE)
		names(popCoarse) = 'offsetExp'
		logPop = log(popCoarse)
		names(logPop) = 'offset'
	} else {
		logPop = NULL
	}
	
	if(!missing(covariates)) {	
		if(!is.list(covariates)) {
			covariates = list(covariates)
			names(covariates) = names(covariates[[1]])[1]
		}
		
		if(verbose) cat("\nrasterizing covariates")
		covariatesFine = geostatsp::stackRasterList(
				covariates,
				rasterFine
		)
		if(verbose) cat("\ndone\n")
		
		if(!missing(border)) {
			covariatesFine = mask(covariatesFine, border)
		}
		if(!missing(population)) {
			covariatesFine = mask(covariatesFine, popFine, maskvalue=0, updatevalue=NA)		
		}
		covariatesCoarse = aggregate(covariatesFine, fact, fun=mean, na.rm=TRUE)
		names(covariatesCoarse) = names(covariates)
	} else {
		covariatesCoarse = NULL
	}
	
	
	zRaster = brick(
			intercept,
			covariatesCoarse,
			logPop
	)
	
	zArray = aperm(as.array(zRaster),c(2,1,3))
	zArray = zArray[,dim(zArray)[2]:1,,drop=FALSE]
	
	zMatrix = matrix(zArray, ncol=dim(zArray)[3])
	theNA = is.na(apply(zMatrix, 1,sum))
	
	zDf = as.data.frame(zMatrix[!theNA,,drop=FALSE])
	zDf$X = 1
	
	zMatrix[theNA,] = 0		
	colnames(zMatrix) = names(zRaster)
	
	anymiss = rep(FALSE, sum(zMatrix[,'intercept']))
	rownames(zDf) = names(anymiss) = as.character(which(zMatrix[,'intercept']>0))
	
	
	attributes(zMatrix) = c(
			attributes(zMatrix),
			polyolay[c('gridobj','fftpoly','mcens','ncens','cellwidth','ext','inclusion')],
			polyolay$gridobj[c('cellInside')],
			list(
					data.frame=zDf,
					M=ncol(zRaster),
					N=nrow(zRaster),
					anymiss=anymiss,
					class=c('lgcpZmat','matrix'),
					pixelOverlay = polyolay$gridobj$pixol,
					missingind=rep(0,ncol(zMatrix))
			)
	)
# haven't added FORM data.frame
	
	
	list(Zmat = zMatrix, events = eventsPPP,
			ext = attributes(zMatrix)$ext,
			cellwidth = attributes(zMatrix)$cellwidth,
			raster = zRaster,
			formula = stats::as.formula(
					paste(
							'X ~ ',
							paste(names(covariates), collapse='+'),
							'+offset(offset)'
					)
			)
	)
}

#' @export
zMatToRaster= function(x) {
	
	Zmat2 = x
	
	gsize = attributes(Zmat2)$cellwidth/2
	
	
	zExtent = extent(
			attributes(Zmat2)$mcens[1]-gsize,
			max(attributes(Zmat2)$mcens)+gsize,
			attributes(Zmat2)$ncens[1]-gsize,
			max(attributes(Zmat2)$ncens)+gsize
	)
	
	rasterMask = raster(
			t(attributes(x)$gridobj$cellInside[
							1:attributes(Zmat2)$M,
							seq(ncol(attributes(x)$gridobj$cellInside),by=-1, 
									len = attributes(Zmat2)$N)
					])	
	)
	extent(rasterMask) = zExtent
	names(rasterMask) = 'mask'
	
	
	zArray2 = array(
			Zmat2, c(
					c(attributes(Zmat2)$M, attributes(Zmat2)$N),#/attributes(Zmat2)$ext,	
					ncol(Zmat2)	
			),
			dimnames=list(NULL, NULL, colnames(Zmat2))
	)
	
	zArray2 = zArray2[,dim(zArray2)[2]:1,]
	zArray2 = aperm(zArray2, c(2,1,3))
	
	zBrick2 = brick(
			zArray2
	)
	extent(zBrick2) = zExtent		
	
	zBrick2 = addLayer(zBrick2, rasterMask)
	
	zBrick2
}


#' @export
pololayToRaster = function(x) {
	polyolay = x
	
	gsize = polyolay$cellwidth/2
	
	rasterMaskExt = raster(
			t(polyolay$gridobj$cellInside[,ncol(polyolay$gridobj$cellInside):1]),
			polyolay$gridobj$mcens[1]-gsize,
			max(polyolay$gridobj$mcens)+gsize,
			polyolay$gridobj$ncens[1]-gsize,
			max(polyolay$gridobj$ncens)+gsize
	)
	zExtent = extent(
			polyolay$gridobj$mcens[1]-gsize-1,
			mean(polyolay$gridobj$mcens),
			polyolay$gridobj$ncens[1]-gsize-1,
			mean(polyolay$gridobj$ncens)
	)		
	
	rasterMask = crop(rasterMaskExt, zExtent)
	
	list(small=rasterMask, ext=rasterMaskExt)
}

