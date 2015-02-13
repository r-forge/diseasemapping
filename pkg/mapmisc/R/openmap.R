osmTiles = function(name) {
	result = c(
			osm = "http://tile.openstreetmap.org",
			"osm-no-labels"="http://a.tiles.wmflabs.org/osm-no-labels/",
      "osm-de"="http://c.tile.openstreetmap.de/tiles/osmde",
			"osm-transport"="http://tile2.opencyclemap.org/transport/",
			"bw-mapnik"="http://b.tiles.wmflabs.org/bw-mapnik2",
			mapquest="http://otile1.mqcdn.com/tiles/1.0.0/osm/",
			"mapquest-sat"="http://otile1.mqcdn.com/tiles/1.0.0/sat",
      "mapquest-labels"='http://otile3.mqcdn.com/tiles/1.0.0/hyb/',
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
#	mapbox="http://a.tiles.mapbox.com/v3/examples.map-vyofok3q/",
	humanitarian="http://a.tile.openstreetmap.fr/hot/",
cartodb='http://c.basemaps.cartocdn.com/light_all/',
'cartodb-dark'='http://c.basemaps.cartocdn.com/dark_all/',
historical='http://www.openhistoricalmap.org/ohm_tiles/'
	)
	
	

  

	# language labels don't appear to be working
	languages = c("en","fr","de", "it","es","ru")
	toadd =	paste("http://a.www.toolserver.org/tiles/osm-labels-", languages,"/", sep="")
	names(toadd) = paste("osm-labels-", languages, sep="")
#	result = c(result, toadd)
	
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
	crs=NA,
  extend=0,
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

	extLL = .extentLL(x,crs, extend)

		
	if(any(abs(as.vector(extLL))>181))
		warning("no CRS supplied and coordinates do not appear to be long-lat")
		
	xlim= c(xmin(extLL), xmax(extLL))
	ylim = c(ymin(extLL), ymax(extLL))
		
	if(missing(zoom)) {
	zoom = 1
	while(nTiles(xlim, ylim, zoom) <= maxTiles & zoom <= 18) {
		zoom = zoom + 1
	}
	zoom = min(c(18,max(c(1, zoom-1))))
	}
	if(verbose) cat("zoom is ", zoom, ", ", nTiles(xlim, ylim, zoom), "tiles\n")

	result = NULL
	for(Dpath in rev(path)) {
		thistile = try(
				getTiles(xlim,ylim, zoom=zoom,
				path=Dpath,
				maxTiles=maxTiles,verbose=verbose),
		silent=TRUE	)

		if(class(thistile)=="try-error"){
			message(paste(Dpath, "not accessible"))
		}	else {
			if(length(names(thistile))) {
				theprefix=strsplit(names(thistile), "([rR]ed|[gG]reen|[bB]lue)$",fixed=FALSE)[[1]]
				names(thistile) = gsub(theprefix, paste(pathOrig[Dpath], "",sep=""), 
					names(thistile),fixed=TRUE)		
			}

		ctable = NULL
		if(!is.null(thistile)) {
			if(nlayers(thistile)==1)
				ctable = thistile@legend@colortable
      result =  stack(thistile, result)	
		}
		
		if(length(ctable))
				result[[1]]@legend@colortable = ctable
		
		} # end not try-error
	} # end loop through path	

	
	if(is.null(result)) {
		result = raster(extLL,1,1,crs=crsLL)
		values(result) = NA
	} 

	crsOut=crs
	if(is.na(crsOut))
		crsOut = projection(x)
	
	if(!is.na(crsOut)  ){

		resultProj = projectRaster(result, crs=crsOut, method="ngb")

    
	} else {
		resultProj = result
	}

		resultProj = stack(resultProj)

		#	resultProj@legend@colortable = result@legend@colortable


	for(D in names(resultProj)) {
			resultProj[[D]]@legend@colortable =
					result[[D]]@legend@colortable
      # set NA's to transparent
      resultProj[[D]]@legend@colortable[1] =
          '#FFFFFF00'
	}

	if(nlayers(resultProj)==1) 
		resultProj = resultProj[[1]]
		

	
	resultProj
}

