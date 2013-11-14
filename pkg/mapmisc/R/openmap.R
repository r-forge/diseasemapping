osmTiles = function(name) {
	result = c(
			osm = "http://tile.openstreetmap.org",
			"osm-no-labels"="http://www.toolserver.org/tiles/osm-no-labels/",
			"osm-transport"="http://tile2.opencyclemap.org/transport/",
			"bw-mapnik"="www.toolserver.org/tiles/bw-mapnik",
			mapquest="http://otile1.mqcdn.com/tiles/1.0.0/osm/",
			"mapquest-sat"="http://mtile02.mqcdn.com/tiles/1.0.0/vy/sat/",
#			hill="http://www.toolserver.org/~cmarqu/hill/",
			landscape="http://tile.opencyclemap.org/landscape/",
#		"osm-retina"="http://tile.geofabrik.de/osm_retina/",
		"opentopomap" = "http://opentopomap.org/tiles/",
#		"osm2world"="http://tiles.osm2world.org/osm/pngtiles/n/",
#		bvg="http://mobil.bvg.de/tiles/",
#	esri="http://server.arcgisonline.com/ArcGIS/rest/services/World_Street_Map/MapServer/tile/",
#	"esri-sat"="http://server.arcgisonline.com/ArcGIS/rest/services/World_Imagery/MapServer/tile/",
#	landshaded="http://tiles.openpistemap.org/landshaded/",
	"maptoolkit"="http://tile2.maptoolkit.net/terrain/",
#	skobbler="http://tiles.skobbler.net/osm_tiles2/",	
	waze="http://tilesworld.waze.com/tiles/",
#	eu="http://alpha.map1.eu/tiles/"#,
#	mapbox="http://a.tiles.mapbox.com/v3/mapbox.natural-earth-hypso-bathy/"
	humanitarian="http://a.tile.openstreetmap.fr/hot/"
	)
	
	
	toolserver = c("parking-bw", "osm-locale-de","bw-noicons")
	toadd =	paste("www.toolserver.org/tiles/", toolserver,"/", sep="")
	names(toadd) = toolserver
	result = c(result, toadd)
	
	
	languages = c("en","fr","de", "it","es","ru")
	toadd =	paste("www.toolserver.org/tiles/osm-labels-", languages,"/", sep="")
	names(toadd) = paste("osm-labels-", languages, sep="")
	result = c(result, toadd)
	
	stamen = c("toner","watercolor")#,"terrain","terrain-background")
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
	pathOrig = path
	pathisname = path %in% names(alltiles)
	path[pathisname] = alltiles[path[pathisname]]
	
	if(length(grep("/$", path, invert=TRUE)))
		path[ grep("/$", path, invert=TRUE)] =
				paste(path[ grep("/$", path, invert=TRUE)], "/", sep="")

	if(length(grep("^http[s]*://", path, invert=TRUE)))
		path[ grep("^http[s]*://", path, invert=TRUE)] = 
				paste("http://", 
						path[ grep("^http[s]*://", path, invert=TRUE)], sep="")
	names(pathOrig) = path
 	
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

	
	# try to get an extent
	
	
	xextent = try(extent(x), silent=TRUE)
	if(class(xextent)=="try-error") {
		# no extent. if it's a vector
		# of length 2, assumit it's a point
		if(is.vector(x)) {
			if(length(x)==2) {
				xextent = extent(x[c(1,1,2,2)])
			} else {
				 warning("if x is a vector extent(x[c(1,1,2,2)]")
			}
		} else {
			warning("can't get an extent from x")
		}
		
	}
	
	# do this because bbox(mybbox) != mybbox
	# but bbox(extent(mybbox) = mybbox
	x = bbox(xextent)
	x2 = x = SpatialPoints(t(x), proj4string=crsIn)
	x = bbox(spTransform(x, CRS("+init=epsg:4326")))
	
	xlim= x[1,]
	ylim = x[2,]
		
	if(missing(zoom)) {
	zoom = 1
	while(nTiles(xlim, ylim, zoom) <= maxTiles & zoom <= 18) {
		zoom = zoom + 1
	}
	zoom = min(c(18,max(c(1, zoom-1))))
	}
	if(verbose) cat("zoom is ", zoom, ", ", nTiles(xlim, ylim, zoom), "tiles\n")

	result = NULL
	

	
	for(Dpath in path) {
		thistile = try(
				getTiles(xlim,ylim, zoom=zoom,
				path=Dpath,
				maxTiles=maxTiles,verbose=verbose),
		silent=TRUE	)

		if(class(thistile)=="try-error"){
			message(paste(Dpath, "not accessible"))
		}	else {
		theprefix=strsplit(names(thistile), "([rR]ed|[gG]reen|[bB]lue)$",fixed=FALSE)[[1]]
		names(thistile) = gsub(theprefix, paste(pathOrig[Dpath], "",sep=""), 
				names(thistile),fixed=TRUE)		
		if(length(result)) {			
			result =  stack(result, thistile)	
		} else {
			result = thistile
		}
	} # end not try-error
	}	

	if(is.null(result)) return(thistile)

	if(all.equal(CRS(proj4string(result)), crs)==TRUE) {
		resultProj = result
	} else {
		resultProj = projectRaster(result, crs=crs, method="ngb")
		# now trim to original bounding box
		pointsOld = SpatialPoints(t(bbox(result)), 
				proj4string=CRS(proj4string(result)))
		pointsNew = spTransform(pointsOld, 
				CRS(proj4string(resultProj)))
		resultProj = crop(resultProj, pointsNew)
	}

	resultProj = stack(resultProj)

	for(D in 1:nlayers(resultProj)) {
		thelen = length(result[[D]]@legend@colortable)
		if(thelen) {
			thevalues = values(resultProj[[D]])
			thevalues[is.na(thevalues)] = thelen
			values(resultProj[[D]]) = thevalues
			
			resultProj[[D]]@legend@colortable =
				c(result[[D]]@legend@colortable, NA)
		}
	}
	
	names(resultProj) = names(result)
	
	if(nlayers(resultProj)==1) {
		resultProj = resultProj[[1]]
	}
	
	
	resultProj
}

