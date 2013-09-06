
openmap = function(upperLeft, lowerRight=NULL, zoom = NULL,
		type = c("osm", "osm-bw", "maptoolkit-topo", "waze", "mapquest", "mapquest-aerial", "bing", "stamen-toner", "stamen-terrain", "stamen-watercolor", "osm-german", "osm-wanderreitkarte", "mapbox", "esri", "esri-topo", "nps", "apple-iphoto", "skobbler", "cloudmade-<id>", "hillshade", "opencyclemap", "osm-transport", "osm-public-transport", "osm-bbike", "osm-bbike-german"),
		minNumTiles = 9L, mergeTiles = TRUE) {

	if(is.vector(upperLeft) ) {
		theproj = NULL
		
	} else {
		
		# do this because bbox(mybbox) != mybbox
		# but bbox(extent(mybbox) = mybbox
		theproj = try(proj4string(upperLeft),silent=TRUE)
		upperLeft = bbox(extent(upperLeft))
		
		if(class(theproj)=="try-error")
			theproj = NULL
 
		
		if(!is.null(theproj)) {
			# transform to long-lat
			upperLeft = SpatialPoints(t(upperLeft),
					proj4string=CRS(theproj))
 
			upperLeft = bbox(spTransform(
							upperLeft, CRS("+proj=longlat")
							))			
		}
		
		lowerRight = c(upperLeft[2,"min"],upperLeft[1,"max"])
		upperLeft = c(upperLeft[2,"max"],upperLeft[1,"min"])
		
	}
 
	
	result = OpenStreetMap::openmap(
			upperLeft,lowerRight, 
			zoom = NULL,
			type = c("osm", "osm-bw", "maptoolkit-topo", "waze", "mapquest", "mapquest-aerial", "bing", "stamen-toner", "stamen-terrain", "stamen-watercolor", "osm-german", "osm-wanderreitkarte", "mapbox", "esri", "esri-topo", "nps", "apple-iphoto", "skobbler", "cloudmade-<id>", "hillshade", "opencyclemap", "osm-transport", "osm-public-transport", "osm-bbike", "osm-bbike-german"),
			minNumTiles = 9L, mergeTiles = TRUE)
	
	if(!is.null(theproj))
		result = OpenStreetMap::openproj(
				result, projection=theproj)
	result
}