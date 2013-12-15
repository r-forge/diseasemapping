

GaussRF = function(x,param=c(variance=1, range=1, shape=1), ...) {
	UseMethod("GaussRF")
	
}

GaussRF.Raster = function(x,param=c(variance=1, range=1, shape=1), ...){


	res = raster(
			RandomFields::RFsimulate(
					modelRandomFields(param),
					as(x, "GridTopology"),
					...
					)
			)
	projection(res) = projection(x)			
	res		
}

GaussRF.SpatialPoints= GaussRF.SpatialPointsDataFrame = 
		function(x,param=c(variance=1, range=1, shape=1), ...){
	
	res = GaussRF.default(
			x=coordinates(x),
			param=param, ...)

	res
}



GaussRF.default = function(x,param=c(variance=1, range=1, shape=1),  ...){
	
	
 	theArgs = list(...)
	theArgs$x = x
	
	if(!any(names(theArgs)=="model")){
		# param is for geostatsp, not RandomFields 

		model  = modelRandomFields(param)
		
		theArgs$model = model	
	} else {
		if(!is.list(theArgs$model)) # if model is a list, it's the extended model definition which doesnt need the param argument
			theArgs$param = param
	}
	
 
	
	result = do.call( RandomFields::RFsimulate, theArgs)

	# some things break (such as spplot) if I add this as an attribute
	#attributes(result)$param = param
	
	res =data.frame(result)
	res[,grep("^variable[[:digit:]]+$", colnames(res))]
	
}


