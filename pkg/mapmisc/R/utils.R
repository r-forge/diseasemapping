# crsLL = CRS("+epsg:4326")
crsLL = CRS("+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0") 

#  crsMerc =CRS("+init=epsg:3857") # mercator projection
crsMercWeb = CRS("+proj=merc +a=6378137 +b=6378137 +lat_ts=0.0 +lon_0=0.0 +x_0=0.0 +y_0=0 +k=1.0 +ellps=WGS84 +units=m +nadgrids=@null +no_defs")
crsMercSimple = CRS("+proj=merc +lat_ts=0.0 +lon_0=0.0 +x_0=0.0 +y_0=0 +k=1.0 +units=m +no_defs +ellps=sphere")

crsMerc = crsMercWeb


openmapExtentLL = extent(-180, 180,-85.0511,85.0511)

if(FALSE) {
  worldLimsLL = SpatialPoints(
      t(bbox(openmapExtentLL)), 
      proj4string=crsLL)
  as.vector(extent(spTransform(worldLimsLL, crsMerc)))
} 

openmapExtentMerc = extent(-20037508,  20037508, -19994838,  19994838)


.getExtentLL = function(x, crs=NA, extend=0) {

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
  
  # force extent to be in allowable bounds
  ymin(result) = pmax(ymin(result), ymin(openmapExtentLL))
  ymax(result) = pmin(ymax(result), ymax(openmapExtentLL))
  
	result
}



.getExtentMerc = function(x, crs=NA, extend=0) {
  
  # find the CRS's
  crsUse = projection(x)
  
  if(all(is.na(crs))) {
    crs=crsLL
  }

  if(all(is.na(crsUse))) {
    crsUse=crs
  }
  
  if(identical(crsUse ,"NA"))
    crsUse = crs
  
  # if x is numeric, transform to extent  
  eps = 1
  if(is.numeric(x))
    if(length(x)==2)
      x = extent(x[1]-eps, xmax=x[1]+eps, ymin=x[2]-eps, ymax=x[2]+eps)
  
  if(requireNamespace('rgdal', quietly=TRUE)) {
    x = as(
        extend(extent(x), extend),
        'SpatialPoints'
    )
    crs(x) = crsUse
    result = extent(spTransform(x, crsMerc))
  } else {
    # no rgdal, try to use raster, doesn't always work
    x = raster(extent(x), nrows=100,ncols=100, crs=crsUse)
    x = raster::extend(x, extend(extent(x), extend))
    result = extent(projectExtent(x, crsMerc))
    
  }
  
  result
}