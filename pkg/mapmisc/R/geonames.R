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

	result = geonames::GNcities(north=north,east=east,
			south=south,west=west,lang,maxRows)

	result = SpatialPointsDataFrame(result[,c("lng","lat")], data=result, 
			proj4string=crsLL)

	if( !identical(projection(theproj), "NA") & ! identical(projection(theproj), NA)) {
		havegdal = require(rgdal, quietly=TRUE )
		if(havegdal)
			result = spTransform(result, CRSobj=CRS(theproj))
	}
		
	result
}

GNsearch = function(...) {
	
	result=geonames::GNsearch(...)
	
	if(all(c("lat","lng") %in% names(result))){
		coords = result[,c("lng","lat"),drop=FALSE]
		mode(coords) = 'numeric'
		result = SpatialPointsDataFrame(
				coords,
				 data=result, 
				proj4string=crsLL)
	}
	result
}