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
			x = extent(x[1], xmax=x[1]+eps, ymin=x[2], ymax=x[2]+eps)

	

	if(requireNamespace('rgdal', quietly=TRUE)) {
    
    # if Spatial object, spTransform it
    
    if(length(grep("^Spatial", class(x)))){
      
      result = extent(
          spTransform(x, crsLL)
          )
      
    } else {
    
      # otherwise, project the extent
    
      x = raster(extent(x), nrows=100,ncols=100, crs=crsUse)
		  result = extent(projectExtent(x, crsLL))
      
    }
    
    # convert extend to ll scale
    if(extend!=0) {
      somePoints = SpatialPoints(
          t(bbox(result)[,c(1,1)]),
          proj4string=crsLL
        )
      somePoints = spTransform(somePoints, CRS(crsUse))
      somePoints@coords[2,] = somePoints@coords[2,] + extend 
      somePoints = spTransform(somePoints, crsLL)
      toextend = as.numeric(apply(somePoints@coords, 2, diff))/2
      
      # make sure not to crop more than possible
      toextend = pmax(toextend,
          - apply(bbox(result),1,diff)/2
          )
      
      result = raster::extend(
          x=result, 
          toextend)
    }
	} else {
    # can't reproject
    
		result = extent(x)
	}	
	result
}