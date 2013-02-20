stackRasterList = function(x, template=x[[1]],method='ngb') {

	if(is.list(x)) {
		if(is.null(names(x)))
			names(x) = paste("c", seq(1, length(x)),sep="")	
	}
	
	Nlayers = length(names(x))
	
	if(length(method)==Nlayers) {
		if(length(names(x)) & all( names(x)%in%names(method)))
			method = method[names(x)]
	} else {
		method = rep(method, Nlayers)
	}
 
			
	result = template
	template = as(template, "BasicRaster")
	for(D in 1:Nlayers) {
		if(as(x[[D]], 'BasicRaster')==template) {
			# same projectoin, same resolution
			result = addLayer(result, x[[D]])			
		}	 else {
			# same projection, different resolution
			if(result@crs@projargs == x[[D]]@crs@projargs) {
				result = addLayer(result,	
						raster::resample(x[[D]], result[[1]])#,method=method[D])
				)
			} else {
				# different resolution
			result = addLayer(result,
					projectRaster(x[[D]], template,method=method[D])
			)
		}
		}
	}
	if(Nlayers == (dim(result)[3]-1) )
		result = result[[-1]]
	names(result) = names(x)
	result
}