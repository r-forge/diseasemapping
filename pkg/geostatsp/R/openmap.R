
openmap = function(x, zoom, 
	path="http://tile.openstreetmap.org/",
	maxTiles = 9,
	crs, cacheDir=tempdir(),
	verbose=FALSE) {

	if(!length(grep("/$", path)))
		path = paste(path, "/", sep="")

	if(!length(grep("^http[s]*://", path)))
		path = paste("http://", path, sep="")
	
	crsIn = try(proj4string(x),silent=TRUE)
	if(missing(crs)) {
		if(class(crsIn)=="try-error"){
			crs =crsIn = CRS("+init=epsg:4326")
		}else { # crs missing but we have crsIn
			crs = crsIn
		}
	} else { # have a crs
		if(class(crsIn)=="try-error"){ # but no crsIn
			crsIn = crs
		}
	}
	
	if(is.character(crs))
		crs = CRS(crs)
	if(is.character(crsIn))
		crsIn = CRS(crsIn)

	
	# do this because bbox(mybbox) != mybbox
	# but bbox(extent(mybbox) = mybbox
	x = bbox(extent(x))
	x2 = x = SpatialPoints(t(x), crsIn)
	x = bbox(spTransform(x, CRS("+init=epsg:4326")))
	
	xlim= x[1,]
	ylim = x[2,]
		
	if(missing(zoom)) {
	zoom = 1
	while(nTiles(xlim, ylim, zoom) <= maxTiles) {
		zoom = zoom + 1
	}
	zoom = max(c(1, zoom-1))
	}
	if(verbose) cat("zoom is ", zoom, ", ", nTiles(xlim, ylim, zoom), "tiles\n")

	
	result =  getTiles(xlim,ylim, zoom=zoom,
				path=path[1],
				maxTiles=maxTiles,verbose=verbose)
			
	if(all.equal(CRS(proj4string(result)), crs)==TRUE) {
		resultProj = result
	} else {
		resultProj = projectRaster(result, crs=crs, method="ngb")
	# for some reason a  bunch of NA's around the edges
	#	extras = (dim(resultProj) - dim(result))[1:2]/2
	#	if(any(extras > 0.5)) {
	#		newextent = extend(extent(resultProj), -extras*res(resultProj))
	#		resultProj = crop(resultProj, newextent)
	#	}
		resultProj@legend@colortable = result@legend@colortable		
	}
	
	resultProj
}

