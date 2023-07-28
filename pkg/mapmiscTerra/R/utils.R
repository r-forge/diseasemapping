

.getRasterMerc = function(zoom) {
  N = 2^zoom 
  terra::rast(ext(extentMerc), nrows = N, ncols=N, crs=crsMerc)
}


.getExtent = function(x, crs=NA, extend=0, crsOut = crsMerc) {
  
	if(is.null(x))
		return(NULL)
	
  # find the CRS's
  crsUse = terra::crs(x)
  
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

		# if long-lat coordinates and buffer is large, assume buffer is metres, do in mercator coords
		if(terra::is.lonlat(crsUse) & any(extend > 80)) {

    	testPointsOut = terra::vect(
    		matrix(as.vector(terra::ext(x)), ncol=2),
    		crs = crsUse
    	)
          testPointsOutMerc = project(testPointsOut, crsMerc)
          outExtentMerc = terra::extend(terra::ext(testPointsOutMerc), extend)
          outPointsMerc = vect(matrix(as.vector(outExtentMerc), ncol=2), crsMerc)
          x = project(outPointsMerc, crsUse)


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

