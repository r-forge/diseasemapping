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
	
	gsize = diff(bRes$ncens[1:2])/2
	
	
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
	
	smallZ = array(
			bRes$Z,
			c(bRes$N, bRes$ext, bRes$M, bRes$ext, ncol(bRes$Z))
	)
	smallZ = drop(smallZ[,1,,1,])
	smallZ = matrix(smallZ, ncol=ncol(bRes$Z))
	smallZ[as.vector(bRes$cellInside)==0, ] = NA
	
	fixed = tcrossprod(smallZ, bRes$betarec)
	fixed = array(fixed, 
			c(bRes$M, bRes$N, nrow(bRes$betarec)))
	fixed = fixed[,dim(fixed)[2]:1,]
	fixed = aperm(fixed, c(2,1,3))
	fixed = brick(fixed)
	extent(fixed) = extent(random)
	fixed = setMinMax(fixed)
	
	fixed = writeRaster(fixed, 
			filename = fixedFile,
  		overwrite = file.exists(fixedFile)
	)
	
	intensity = exp(fixed + random)
	intensity = setMinMax(intensity)
	
	intensity = writeRaster(intensity,
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
	intensitySummary = calc(intensity, summaryFunction)
	names(intensitySummary) = gsub("^xx", 'intensity', names(intensitySummary))
	
	summaryStack = stack(
			randomSummary, intensitySummary
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
			intensity = intensity,
			random = random,
			fixed = fixed,
			mask=rasterMask,
			beta = bRes$betarec,
			sd = exp(bRes$etarec[,1]),
			range = 2*exp(bRes$etarec[,2])
	)
	samples$beta[,1] =
			samples$beta[,1] - samples$sd^2/2
	
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
					y = dlnorm(xSeq, 
							bRes$priors$etaprior$mean[1], 
							bRes$priors$etaprior$variance[1,1]
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
					y = dlnorm(xSeq, 
							bRes$priors$etaprior$mean[2], 
							bRes$priors$etaprior$variance[2,2]
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

