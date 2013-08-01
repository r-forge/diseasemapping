GaussRF = function(x,param, ...) {
	UseMethod("GaussRF")
	
}

GaussRF.Raster = function(x,param, ...){

	xseq = seq(x@extent@xmin, x@extent@xmax, len=x@ncols)
	yseq = seq(x@extent@ymin, x@extent@ymax, len=x@nrows)

	thediffx = diff(range(diff(xseq)))
	if(thediffx > 0) {
		thediffx = signif(diff(xseq[1:2]), 10^-ceiling(log10(thediffx)+1))
		xseq = seq(x@extent@xmin, by=thediffx, len=x@ncols )
	}

	thediffy = diff(range(diff(yseq)))
	if(thediffy > 0) {
		thediffy = signif(diff(yseq[1:2]), 10^-ceiling(log10(thediffy)+1))
		yseq = seq(x@extent@ymin, by=thediffy, len=x@ncols )
	}
	
	
	
	res = GaussRF(x=xseq, param=param,
			y=yseq,	grid=TRUE, 
		...
			)

			
	resRast = raster(t(res[,seq(dim(res)[2], 1)]),
			x@extent@xmin, x@extent@xmax,
			x@extent@ymin, x@extent@ymax, crs=x@crs)
	
	return(resRast)
}

GaussRF.SpatialPointsDataFrame = function(x,param, ...){
	
	x=coordinates(x)
	NextMethod("GaussRF")
}

GaussRF.SpatialPoints= function(x,param, ...){
	
	x=coordinates(x)
	NextMethod("GaussRF")
}

GaussRF.default = function(x,param,  ...){
	
	
 	theArgs = list(...)
	theArgs$x = x
	
	if(!any(names(theArgs)=="model")){
		# param is for geostatsp, not RandomFields 

		requiredParams = c("variance","range","maternRoughness")
		if(!all(requiredParams %in% names(param)))
			warning("param has names", paste(names(param),collapse=","), 
					" must have ", paste(requiredParams, collapse=","))

		param["scale"] = 2*param["range"] 
		
		if(any(names(param)=="aniso.ratio")){
			# geometric anisotropy
			if(any(names(param)=="aniso.angle.degrees") & 
					!any(names(param)=="aniso.angle.radians") ) {
				param["aniso.angle.radians"] = param["aniso.angle.degrees"]*2*pi/360				
		}
			model=list("$", var=param["variance"],   
					A=anisoMatrix(-param["aniso.angle.radians"],
							param["scale"]*c(param["aniso.ratio"],1)
				),
				list("matern", nu=param["maternRoughness"]))	
		} else {
		
			model=list("$", var=param["variance"],   
					s=param["scale"],
					list("matern", nu=param["maternRoughness"]))	
			
		}	
		theArgs$model = model	
	} else {
		theArgs$param = param
	}
	
	
	do.call( RandomFields::GaussRF, theArgs) 
	#(x, model=model, param=param,	...)
	
	
}


