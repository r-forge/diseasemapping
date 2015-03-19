# crsLL = CRS("+epsg:4326")
crsLL = CRS("+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0") 

.extentLL = function(x, crs=NA, extend=0) {

  # find the CRS's
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
  
  # if x is numeric, transform to extent  
  eps = 1e-5
	if(is.numeric(x))
		if(length(x)==2)
			x = extent(x[1]-eps, xmax=x[1]+eps, ymin=x[2]-eps, ymax=x[2]+eps)


  if(requireNamespace('rgdal', quietly=TRUE)) {
  x = as(
      extend(extent(x), extend),
      'SpatialPoints'
      )
  crs(x) = crsUse
  result = extent(spTransform(x, crsLL))
} else {
  # no rgdal, try to use raster, doesn't always work
    x = raster(extent(x), nrows=100,ncols=100, crs=crsUse)
    x = raster::extend(x, extend(extent(x), extend))
    result = extent(projectExtent(x, crsLL))
    
  }
  
	result
}