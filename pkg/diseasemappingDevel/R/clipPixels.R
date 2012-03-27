clipPixels = function(x, buffer=NULL, poly=NULL) {
	
	if(!is.null(buffer)) {
		therange =  apply(coordinates(x), 2, range)
		therange[1,] = therange[1,]+  buffer
		therange[2,] = therange[2,]-  buffer
		
		x = x[   
				coordinates(x)[,1] > therange[1,1] &
        				coordinates(x)[,1] < therange[2,1] &
						coordinates(x)[,2] > therange[1,2] &
						coordinates(x)[,2] < therange[2,2] ,
		]
		
		
	}
	
	if(!is.null(poly)) {
		inWater = is.na(overlay(x, spTransform(poly, CRS(proj4string(x)))))
		x@data[inWater, ] = NA
	}
	
	x
	
}