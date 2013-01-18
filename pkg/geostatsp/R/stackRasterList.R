stackRasterList = function(x, template=x[[1]],method='ngb') {

	if(length(method)==length(x)) {
		if(all( names(x)%in%names(method)))
			method = method[names(x)]
	} else {
		method = rep(method, length(x))
	}
	names(method) = names(x)
			
	result = template
	template = as(template, "BasicRaster")
	for(D in names(x)) {
		if(as(x[[D]], 'BasicRaster')==template) {
			result = addLayer(result, x[[D]])			
		}	 else {
			result = addLayer(result,
					projectRaster(x[[D]], template,method=method[[D]])
			)
		}
	}
	if(length(x) == (length(names(result))+1) )
		result = result[[-1]]
	names(result) = names(x)
	result
}