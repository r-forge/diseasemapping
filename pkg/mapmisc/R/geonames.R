GNcities = function(north, east, south, west, lang = "en", maxRows = 10) {
	
	fourCoords=FALSE
	if(is.numeric(north))
		if(length(north)== 1)
			fourCoords = TRUE

	theproj = projection(north)
	if(!fourCoords) {
		extLL = .extentLL(north)
		
		east = xmax(extLL)
		west = xmin(extLL)
		south = ymin(extLL)
		north= ymax(extLL)
		
	}

	if (requireNamespace("geonames", quietly = TRUE)) { 
		
	result = geonames::GNcities(north=north,east=east,
			south=south,west=west,lang,maxRows)

	result = SpatialPointsDataFrame(cbind(
					as.numeric(result[,'lng']),
					as.numeric(result[,'lat'])
					), data=result, 
			proj4string=crsLL)
} else {
	warning("install the geonames package to use GNcities")
	result = NULL
}

	if( !identical(projection(theproj), "NA") & ! identical(projection(theproj), NA)) {
		if(requireNamespace('rgdal', quietly=TRUE ))
			result = spTransform(result, CRSobj=CRS(theproj))
	}
		
	result
}

GNsearch = function(...) {
	
	
	if(requireNamespace("geonames", quietly = TRUE)) { 
	result=geonames::GNsearch(...)
	
	if(all(c("lat","lng") %in% names(result))){
		coords = as.matrix(result[,c("lng","lat"),drop=FALSE])
		mode(coords) = 'numeric'

		result$population = as.numeric(result$population)
		
		result = SpatialPointsDataFrame(
				coords,
				 data=result, 
				proj4string=crsLL)
	}
} else {
	warning("install the geonames package to use GNsearch")
	result = NULL
	
}
	result
}