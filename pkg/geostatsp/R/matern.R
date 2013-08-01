matern = function(x, y=x, range, maternRoughness, variance=1, 
		aniso.ratio=1, aniso.angle.degrees=0, 
		aniso.angle.radians=aniso.angle.degrees*2*pi/360) {
	UseMethod("matern")
	
}

matern.dist = function( x, y=x, range, maternRoughness, variance=1, 
		aniso.ratio=1, aniso.angle.degrees=0, 
		aniso.angle.radians=aniso.angle.degrees*2*pi/360) {
	x = as.matrix(x)
	NextMethod("matern")
}

matern.SpatialPointsDataFrame = function( x, y=x, range, maternRoughness, variance=1, 
		aniso.ratio=1, aniso.angle.degrees=0, 
		aniso.angle.radians=aniso.angle.degrees*2*pi/360) {
	x = SpatialPoints(x)
	NextMethod("matern")
}


matern.Raster = function( x, y=x, range, maternRoughness, variance=1, 
		aniso.ratio=1, aniso.angle.degrees=0, 
		aniso.angle.radians=aniso.angle.degrees*2*pi/360)
 {

	 if(length(grep("^Raster", class(y)))){
		 y = SpatialPoints(as.data.frame(y, xy=TRUE)[,c("x","y")])
	 }
	 
	resRast = x
	x= SpatialPoints(as.data.frame(x, xy=TRUE)[,c("x","y")])

	# for some reason
	# NextMethod("matern")
	# doesn't work, it calls matern.default
	result= matern(x, y, range, maternRoughness, variance, 
		aniso.ratio, aniso.angle.degrees, 
		aniso.angle.radians)

	if(length(result)==ncell(resRast)) {
			values(resRast) = result 
			result=resRast
	} 
	result
}



matern.SpatialPoints = function( x, y=x, range, maternRoughness, variance=1, 
		aniso.ratio=1, aniso.angle.degrees=0, 
		aniso.angle.radians=aniso.angle.degrees*2*pi/360
		){
		if(length(grep("SpatialPoints", class(y)))) {
			y = y@coords[,1] + 1i*y@coords[,2]  
		}
		if(length(y)==2 & !is.complex(y)){
			y = y[1] + 1i*y[2]
		}
		x = x@coords[,1] + 1i*x@coords[,2]
		
		thedist = outer(x, y, FUN="-")
		if(aniso.ratio != 1){ # geometric anisotropy
			thedist = thedist * exp(1i*aniso.angle.radians)
			thedist =Re(thedist) +  (1i/aniso.ratio )*Im(thedist)
		}
		x = Mod(thedist)
		
		# for some reason
		# NextMethod("matern")
		# doesn't work, it calls matern.default
		result= matern(x, y, range, maternRoughness, variance, 
				aniso.ratio, aniso.angle.degrees, 
				aniso.angle.radians)
		
		attributes(result)$aniso = c(ratio=aniso.ratio, angle.radians=aniso.angle.radians)
		result
}

matern.default = function( x, y=x, range, maternRoughness, variance=1, 
		aniso.ratio=1, aniso.angle.degrees=0, 
		aniso.angle.radians=aniso.angle.degrees*2*pi/360)
{
	
	
	
	xscale = abs(x)*(sqrt(8*maternRoughness)/ range)
	result = ( variance/(gamma(maternRoughness)* 2^(maternRoughness-1)  ) ) * 
			( xscale^maternRoughness *
				besselK(xscale , maternRoughness) )
	result[x==0] = variance
	
	attributes(result)$model = c(range=range, maternRoughness=maternRoughness,
			variance=variance)
	result
	
}
