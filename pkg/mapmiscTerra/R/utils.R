
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

extentMerc = terra::ext(-20037508.3427892, 20037508.3427892, -20037508.3427893, 20037508.3427892)

.getRasterMerc = function(zoom) {
  N = 2^zoom 
  terra::rast(extentMerc, nrows = N, ncols=N, crs=crsMerc)
}


.getExtent = function(x, crs=NA, extend=0, crsOut = crsMerc) {
  
	if(is.null(x))
		return(NULL)
	
  # find the CRS's
  crsUse = crs(x)
  
#	if(is.logical(crs)) { # it's probably NA,
  if(all(is.na(crs))) {
    crs=crsLL
  }

  if(is.na(crsUse)) crsUse = crs
  
  # if x is numeric, transform to extent  
  tileEps = sqrt(.Machine$double.neg.eps)
  if(is.numeric(x)) {
    if(length(x)==2) {
      x = terra::ext(x[1]-tileEps, x[1]+tileEps, x[2]-tileEps, x[2]+tileEps)
		}
	}

		# if long-lat coordinates and buffer is large, assume buffer is metres 
		if(terra::is.lonlat(crs) & any(extend > 80)) {
    	x = terra::vect(
    		matrix(as.vector(terra::ext(x)), ncol=2),
    		crs = crsUse
    	)


    	x = terra::vect(
	  		geosphere::destPoint(
    			    		matrix(as.vector(terra::ext(x)), ncol=2),
									c(-45,45, 135,-135),
									extend * sqrt(2)
    			),
    		crs = crsLL)


		} else {
				x = vect(
						matrix(as.vector(terra::extend(terra::ext(x), extend)), ncol=2),
						crs = crsUse)

		}

    result =  terra::ext(terra::project(x, crsOut))

  
    result
}


.cropExtent = function(x,y){
	x = as.vector(x)
	res = as.vector(y)
	res[c(1,3)] = pmin(res[c(1,3)], x[c(1,3)])
	res[c(2,4)] = pmax(res[c(2,4)], x[c(2,4)])

  terra::ext(x)
}

