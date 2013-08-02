matern = function(x, y=x, param=c(range=1, variance=1, rough=1)) {
	UseMethod("matern")
	
}

matern.dist = function( x, y=x, param=c(range=1, variance=1, rough=1)) {
	x = as.matrix(x)
	NextMethod("matern")
}

matern.SpatialPointsDataFrame = function( x, y=x, param=c(range=1, variance=1, rough=1)) {
	x = SpatialPoints(x)
	NextMethod("matern")
}


matern.Raster = function( x, y=x, param=c(range=1, variance=1, rough=1))
 {

	 if(length(grep("^Raster", class(y)))){
		 y = SpatialPoints(as.data.frame(y, xy=TRUE)[,c("x","y")])
	 }
	 
	resRast = x
	x= SpatialPoints(as.data.frame(x, xy=TRUE)[,c("x","y")])

	# for some reason
	# NextMethod("matern")
	# doesn't work, it calls matern.default
	result= matern(x, y, param)

	if(length(result)==ncell(resRast)) {
			values(resRast) = result 
			result=resRast
	} 
	result
}



matern.SpatialPoints = function(x, y=x,param=c(range=1, variance=1, rough=1)
		){
		if(length(grep("SpatialPoints", class(y)))) {
			y = y@coords[,1] + 1i*y@coords[,2]  
		}
		if(length(y)==2 & !is.complex(y)){
			y = y[1] + 1i*y[2]
		}
		x = x@coords[,1] + 1i*x@coords[,2]
		
		thedist = outer(x, y, FUN="-")
		
		if(any(names(param)=="aniso.ratio") ) { # geometric anisotropy
			if(any(names(param)=="aniso.angle.degrees") & 
				!any(names(param)=="aniso.angle.radians") ) {
			param["aniso.angle.radians"] = param["aniso.angle.degrees"]*2*pi/360				
			}
			if(!any(names(param)=="aniso.angle.radians") )
				warning("anisotropy angle not supplied")
			
			thedist = thedist * exp(-1i*param["aniso.angle.radians"])
			thedist =Re(thedist) +  (1i/ param["aniso.ratio"] )*Im(thedist)
		}
		x = Mod(thedist)
		
		# for some reason
		# NextMethod("matern")
		# doesn't work, it calls matern.default
		result= matern(x, y, param)
		
		result
}

matern.default = function( x, y=x,param=c(range=1, variance=1, rough=1))
{
	
	names(param) = gsub("^var$", "variance", names(param))
	if(!any(names(param)=="variance") )
		param["variance"]=1
		
	
	xscale = abs(x)*(sqrt(8*param["rough"])/ param["range"])
	result = ( param["variance"]/(gamma(param["rough"])* 2^(param["rough"]-1)  ) ) * 
			( xscale^param["rough"] *
				besselK(xscale , param["rough"]) )
	result[x==0] = param["variance"]
	
	attributes(result)$param = param
	result
	
}
