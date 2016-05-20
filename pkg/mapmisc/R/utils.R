# crsLL = CRS("+init=epsg:4326")
crsLL = CRS("+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0") 

# CRS("+init=epsg:3857") without the nagrids stuff
crsMerc = CRS("+proj=merc +a=6378137 +b=6378137 +lat_ts=0.0 +lon_0=0.0 +x_0=0.0 +y_0=0 +k=1.0 +units=m +no_defs")

#crsMerc =  CRS("+proj=merc +ellps=sphere +units=m")
#crsLlSphere = CRS("+proj=longlat +ellps=sphere")

# extent of the Earth in the spherical mercator projection
if(FALSE) {
# need rgdal for this	
createExtentMerc  = function(){

	latlim = atan(sinh(pi))*360/(2*pi)
	lonlim = 180

	openmapExtentLL = extent(-lonlim, lonlim,-latlim,latlim)

	extentMerc = extent(projectExtent(raster(openmapExtentLL, crs=crsLL), crs=crsMerc))

	extentMerc
}

extentMerc = createExtentMerc()
dput(extentMerc, file='')
}

extentMerc = new("Extent"
    , xmin = -20037508.3427892
    , xmax = 20037508.3427892
    , ymin = -20037508.3427893
    , ymax = 20037508.3427892
)

.getRasterMerc = function(zoom) {
  N = 2^zoom 
  raster(extentMerc, nrows = N, ncols=N, crs=crsMerc)
}


.getExtent = function(x, crs=NA, extend=0, crsOut = crsMerc) {
  
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
  tileEps = sqrt(.Machine$double.neg.eps)
  if(is.numeric(x)) {
    if(length(x)==2) {
      x = extent(x[1]-tileEps, xmax=x[1]+tileEps, ymin=x[2]-tileEps, ymax=x[2]+tileEps)
		}
	}
  if(requireNamespace('rgdal', quietly=TRUE)) {
		# if long-lat coordinates and buffer is large, assume buffer is metres 
		if(raster::isLonLat(crs) & any(extend > 80)) {
    	x = as(
        	extent(x), 
        	'SpatialPoints'
    	)
    	crs(x) = crsUse
			x@coords = geosphere::destPoint(
					x@coords,
					c(-45,45, 135,-135),
					extend * sqrt(2)
					)
		} else {
			x = as(
					extent(raster::extend(extent(x), extend)),
					'SpatialPoints'
			)
			crs(x) = crsUse
		}

    result =  extent(spTransform(x, crsOut))
  } else {
    # no rgdal, try to use raster, doesn't always work
    x = raster(extent(x), nrows=100,ncols=100, crs=crsUse)
    x = raster::extend(extent(x), extend)
    result = extent(projectExtent(x, crsOut))
  }
  
  result
}


.cropExtent = function(x,y){
  xmin(x) = pmax(xmin(x), xmin(y))
  xmax(x) = pmin(xmax(x), xmax(y))
  ymin(x) = pmax(ymin(x), ymin(y))
  ymax(x) = pmin(ymax(x), ymax(y))
  x
}

