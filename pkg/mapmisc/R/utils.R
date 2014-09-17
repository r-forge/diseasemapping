# crsLL = CRS("+epsg:4326")
crsLL = CRS("+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0") 

.extentLL = function(x, crs=NA) {

	eps = 1e-5
	if(is.numeric(x))
		if(length(x)==2)
			x = extent(x[1], xmax=x[1]+eps, ymin=x[2], ymax=x[2]+eps)

	crsUse = projection(x)

	
#	if(is.logical(crs)) { # it's probably NA,
	if(all(is.na(crs))) {
		crs=crsLL
	}
#	if(is.logical(crsUse)) { # it's probably NA,
	if(all(is.na(crsUse))) {
		crsUse=crs
	}
	
	if(identical(crsUse ,"NA"))
		crsUse = crs
	
	x = raster(extent(x), nrows=100,ncols=100, crs=crsUse)

	if(requireNamespace('rgdal', quietly=TRUE)) {
		result = extent(projectExtent(x, crsLL))
	} else {
		result = extent(x)
	}	
	result
}