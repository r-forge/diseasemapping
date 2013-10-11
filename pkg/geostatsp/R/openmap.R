osmTiles = function(name) {
	result = c(
			osm = "http://tile.openstreetmap.org",
			"osm-no-labels"="http://www.toolserver.org/tiles/osm-no-labels/",
			"osm-transport"="http://tile2.opencyclemap.org/transport/",
			"bw-mapnik"="www.toolserver.org/tiles/bw-mapnik",
			mapquest="http://otile1.mqcdn.com/tiles/1.0.0/osm/"
	)
	
	languages = c("en","fr","de", "it")
	toadd =	paste("www.toolserver.org/tiles/osm-labels-", languages,"/", sep="")
	names(toadd) = paste("osm-labels-", languages, sep="")
	result = c(result, toadd)
	
	stamen = c("terrain","toner","watercolor")
	toadd = paste("http://tile.stamen.com/", stamen, "/",sep="")
	names(toadd) = paste("stamen-", stamen, sep="")
	
	result = c(result, toadd)
	
	
	if(!missing(name)) {
		if(name %in% names(result)) {
			result = result[name]
		} else {
			warning("name ", name, " is not a tile path, returning vector of all tile paths")
		}
	}
	
	result
	
}

openmap = function(x, zoom, 
	path="http://tile.openstreetmap.org/",
	maxTiles = 9,
	crs,  
	verbose=FALSE) {

	alltiles = osmTiles()
	pathisname = path %in% names(alltiles)
	path[pathisname] = alltiles[path[pathisname]]
	
	if(length(grep("/$", path, invert=TRUE)))
		path[ grep("/$", path, invert=TRUE)] =
				paste(path[ grep("/$", path, invert=TRUE)], "/", sep="")

	if(length(grep("^http[s]*://", path, invert=TRUE)))
		path[ grep("^http[s]*://", path, invert=TRUE)] = 
				paste("http://", 
						path[ grep("^http[s]*://", path, invert=TRUE)], sep="")
	
 	
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
	zoom = min(c(18,max(c(1, zoom-1))))
	}
	if(verbose) cat("zoom is ", zoom, ", ", nTiles(xlim, ylim, zoom), "tiles\n")

	result = NULL
	for(Dpath in path) {
		thistile = getTiles(xlim,ylim, zoom=zoom,
				path=Dpath,
				maxTiles=maxTiles,verbose=verbose)
		names(thistile) = gsub("^http://|/$", "", Dpath)
		if(length(result)) {
		result =  stack(result, thistile)	
		} else {
			result = thistile
		}
	}	


	if(all.equal(CRS(proj4string(result)), crs)==TRUE) {
		resultProj = result
	} else {
		resultProj = stack(projectRaster(result, crs=crs, method="ngb"))

		
	  
	}
	
	result <<- result
	resultProj <<- resultProj
	
	for(D in 1:nlayers(resultProj)) {
		thelen = length(result[[D]]@legend@colortable)
		resultProj[[D]] = calc(resultProj[[D]], function(qq) {
					qq[is.na(qq)] = thelen
					qq
				})
		
		resultProj[[D]]@legend@colortable =
				c(result[[D]]@legend@colortable, NA)
	}
	
	names(resultProj) = names(result)
	
	if(nlayers(resultProj)==1) {
		resultProj = resultProj[[1]]
	}
	
	
	resultProj
}

