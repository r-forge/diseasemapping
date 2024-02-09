writeRasterMapTiles = function(x, filename, overwrite = TRUE,  ...) {
    if(!length(filename)){
    	return(x)
    }

  theTiles = attributes(x)$tiles
  theOpenmap = attributes(x)$openmap


  if(any(nrow(terra::coltab(x)[[1]])>256) ) {
        # color table too big for png
        # transparency will be dropped
    filename = gsub("png$", "tif", filename)
    x2 = terra::colorize(x, to = 'rgb', alpha=FALSE)
    result = terra::writeRaster(x2, filename = filename,
      datatype = 'INT1U', overwrite = overwrite, ...)
  } else if(terra::has.RGB(x)) {
    	# save rgb layers as tif
      # drop transparency
    	filename = gsub("png$", "tif", filename)
      result = terra::writeRaster(x,  
        datatype = 'INT1U',
        filename=filename, overwrite=overwrite, ...)
	} else if(terra::has.colors(x)) {
    # write as png
    filename = gsub("tif$", "png", filename)
    # make colour table go to 255
    coltabOrig = terra::coltab(x)[[1]]
    origValues =coltabOrig[,1]
    if(any(origValues > 255 | origValues < 0,na.rm=TRUE )) warning("raster values outside range 0 to 255")
    newValues = setdiff(0:255, origValues)
    if(length(newValues)) {
        coltabToAdd = coltabOrig[rep(1, length(newValues)), ]
        coltabToAdd[,1] = as.integer(newValues)         
        coltabToAdd[,-1] = 0L
        newColtab = rbind(coltabOrig, coltabToAdd)
        newColtab = newColtab[order(newColtab[,1]),]
        terra::coltab(x) = newColtab
    }
    result = terra::writeRaster(x, filename, 
      datatype = 'INT1U', overwrite=overwrite, ...)
  } else {
    warning("no colour table or RGB")
    result = terra::writeRaster(x, filename, 
      overwrite=overwrite, ...)
  }


  attributes(result)$tiles = theTiles
  attributes(result)$openmap = theOpenmap
  attributes(result)$source = terra::sources(result)
  result
}


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

