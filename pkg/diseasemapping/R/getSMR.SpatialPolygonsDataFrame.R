getSMR_SpatialPolygonsDataFrame <- function(popdata, model, casedata=NULL, 
		regionCode = intersect(names(popdata), names(casedata))[1],
		regionCodeCases=regionCode, area=FALSE,  ...) {

	
	if (area & !("surfaceArea" %in% names(popdata) ) ) {
		if(length(grep("longlat", popdata@proj4string@projargs)))
			warning("computing areas of polygons in long-lat projection")
		popdata$surfaceArea = sapply(slot(popdata, "polygons"), slot, "area")
	}	
	
popdata@data <- getSMR(popdata@data, model, casedata, regionCode,
		regionCodeCases, area,...)



popdata
}
