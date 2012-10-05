mrfGridData <-
function(points,Napprox=11, 
		covariates=NULL, offset=NULL, buffer=2) {
	
	
	thebb = bbox(points)
		
	Nx = Napprox
	cellSize = (thebb[1,2]-thebb[1,1])/Nx
	Ny = (thebb[2,2]-thebb[2,1])/cellSize
	thebb[2,2]= thebb[2,1]+cellSize*Ny
	
	rasterTemplate = raster(nrows=Ny, ncols=Nx, 
			xmn=thebb[1,1], xmx=thebb[1,2], 
			ymn=thebb[2,1], ymx=thebb[2,2])	
	
	# vector of counts per cell
	count = rasterize(thesim$points, rasterTemplate,
			field=rep(1,length(thesim$points)),
			fun='sum'
	)
	
	count[][is.na(count[])] = 0
	
	result = data.frame(
			count = count[],
			offset= 2*log(cellSize) 
	)
	
	if(!is.null(offset)) {
		result$offset = result$offset +
				resample(offset,
						rasterTemplate)[]
	}
	
	# resample covariates, convert to vector
	if(!is.null(covariates)) {
		
		covariates@layernames[covariates@layernames==''] = 
				paste('cov', 
						seq(1,sum(covariates@layernames=='')),
						sep='')	
		
		for(Dcov in covariates@layernames) {
			result[[Dcov]] =
					resample(covariates[[Dcov]],
							rasterTemplate)[]
		}
	}	
	
	
	# now the vector of indices
	Nx2 = Nx + 2*buffer
	Ny2 = Ny + 2*buffer
	
	xvec = seq(from=buffer+1, by=1, len=Nx)
	xvec = rep(xvec, Ny )
	yvec= seq(from=buffer+1, by=1, len=Ny)
	yvec =  rep(yvec, rep(Ny, Nx))
	
# indices consistent with rasters.  
# cell 1 is lower left corner,
#  then move along columns.
	
	index = xvec + Nx2*(yvec-1)
	
	
	result$cellID = index
	
	rasterTemplateBig = raster(
				nrows=Ny2, ncols=Nx2, 
					xmn=thebb[1,1]-buffer*cellSize, 
					xmx=thebb[1,2]+buffer*cellSize, 
					ymn=thebb[2,1]-buffer*cellSize, 
					ymx=thebb[2,2]+buffer*cellSize,
					crs=CRS(proj4string(points))
			)	
		
	attributes(result)$rasterTemplate = rasterTemplateBig	
	
	attributes(result)$buffer = buffer
	attributes(result)$bbox = thebb
	attributes(result)$nrowForInla = ncol(rasterTemplateBig)
	attributes(result)$ncolForInla = nrow(rasterTemplateBig)
	
	result	
	
}
