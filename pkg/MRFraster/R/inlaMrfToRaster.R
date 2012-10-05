inlaMrfToRaster <-
function(x) {

	rasterTemplate = attributes(x$.args$data)$rasterTemplate
	bbox = attributes(x$.args$data)$bbox
	
	if(is.null(rasterTemplate)){
		warning("no raster template, was inla run with rasterForMrf?")
	}
	
	mrf = names(x$summary.random)[1]
	
	resRaster = brick(rasterTemplate,
			nl=dim(x$summary.random[[mrf]])[2])

	resRaster[] = as.matrix(x$summary.random[[mrf]])
	if(!is.null(bbox))
		resRaster = crop(resRaster,extent(bbox))

	resRaster@layernames = 
			paste("random.",resRaster@layernames,sep="")
	
	if(!is.null(x$summary.fitted.values)) {
	fittedData = x$summary.fitted.values
	
	theFormula = terms(x$.args$formula)
	theoffset = attributes(theFormula)$offset
	if(!is.null(theoffset)) {
		theoffset = as.character(attributes(theFormula)$variables[[theoffset+1]])[2]
		addOffset = grep("mean|quant", colnames(fittedData))	
		toAdd = matrix(0, nrow(fittedData), ncol(fittedData))
		toAdd[,addOffset] = x$.args$data[,theoffset]
		fittedData = fittedData *exp(- toAdd)		
	}
	
	
	fittedData[,mrf] = x$.args$data[,mrf]
	
	allCells = seq(1,
			length(rasterTemplate))
	NotInFitted = allCells[!allCells %in% fittedData[,mrf]]
	
	newMat = matrix(NA,length(NotInFitted), ncol(fittedData))
	colnames(newMat) = colnames(fittedData)
	newMat[,mrf]=NotInFitted
	fittedData = rbind(fittedData,newMat)
	fittedData=fittedData[order(fittedData[,mrf]),]
	
	
	fittedRaster = 	brick(rasterTemplate,
			nl=ncol(fittedData))
	fittedRaster[]= as.matrix(fittedData)
 
	fittedRaster = crop(fittedRaster,resRaster)
	
	fittedRaster@layernames = paste("fitted.",fittedRaster@layernames,sep="")
	resRaster = stack(resRaster, fittedRaster)
	
} 
	
	
	return(resRaster)
	
}
