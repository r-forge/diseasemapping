osmTiles = function(name) {
	result = c(
			osm = "http://tile.openstreetmap.org",
			'osm-admin' = 'http://korona.geog.uni-heidelberg.de/tiles/adminb',
			'osm-roads-grey' = 'http://korona.geog.uni-heidelberg.de/tiles/roadsg',
			'osm-roads' = 'http://korona.geog.uni-heidelberg.de/tiles/roads',
			'osm-semitransparent' = 'http://korona.geog.uni-heidelberg.de/tiles/hybrid',
			"osm-no-labels"="http://tiles.wmflabs.org/osm-no-labels/",
      "osm-de"="http://c.tile.openstreetmap.de/tiles/osmde",
			"osm-transport"="http://tile2.opencyclemap.org/transport/",
			"bw-mapnik"="http://tiles.wmflabs.org/bw-mapnik",
			mapquest="http://otile1.mqcdn.com/tiles/1.0.0/osm/",
			"mapquest-sat"="http://otile1.mqcdn.com/tiles/1.0.0/sat",
      "mapquest-labels"='http://otile3.mqcdn.com/tiles/1.0.0/hyb/',
      'osm-cyclemap' = 'http://a.tile.opencyclemap.org/cycle/',
      'osm-seamap' = 'http://tiles.openseamap.org/seamark/',
      'osm-fr' = 'http://a.tile.openstreetmap.fr/osmfr/',
			'landscape'="http://tile.opencyclemap.org/landscape/",
			'hyda' = 'http://c.tile.openstreetmap.se/hydda/full/',
			'hyda-base' = 'http://c.tile.openstreetmap.se/hydda/base/',
			'hyda-roads' = 'http://c.tile.openstreetmap.se/hydda/roads_and_labels/',
  		"opentopomap" = "http://a.tile.opentopomap.org/tiles/",
 "maptoolkit"="http://tile2.maptoolkit.net/terrain/",
	waze="http://worldtiles3.waze.com/tiles/",
# 'waze-us'='https://livemap-tiles2.waze.com/tiles/',
	humanitarian="http://a.tile.openstreetmap.fr/hot/",
  cartodb='http://c.basemaps.cartocdn.com/light_all/',
 'cartodb-dark'='http://c.basemaps.cartocdn.com/dark_all/',
  historical='http://www.openhistoricalmap.org/ohm_tiles/',
  nrcan = 
		'http://geoappext.nrcan.gc.ca/arcgis/rest/services/BaseMaps/CBMT_CBCT_GEOM_3857/MapServer/tile/',
 'nrcan-text' = 
		 'http://geoappext.nrcan.gc.ca/arcgis/rest/services/BaseMaps/CBMT_TXT_3857/MapServer/tile/',
 'nrcan-text-fr' = 
		 'http://geoappext.nrcan.gc.ca/arcgis/rest/services/BaseMaps/CBCT_TXT_3857/MapServer/tile/'
)
#	skobbler="http://tiles3.skobbler.net/osm_tiles2/",	
#  'sputnik' = 'http://a.tiles.maps.sputnik.ru/tiles/kmt2/',
#		"osm2world"="http://tiles.osm2world.org/osm/pngtiles/n/",
#		bvg="http://mobil.bvg.de/tiles/",
#	landshaded="http://tiles.openpistemap.org/landshaded/",
#		"osm-retina"="http://tile.geofabrik.de/osm_retina/",
#      'osm-rail' = 'http://a.tiles.openrailwaymap.org/standard/',
# rail is 512 insstead of 256 tiles
#			hill="http://www.toolserver.org/~cmarqu/hill/",
#	eu="http://alpha.map1.eu/tiles/",
# 'esri' = 'http://services.arcgisonline.com/ArcGIS/rest/services/World_Street_Map/MapServer/tile/',
# 'esri-grey' = 'http://services.arcgisonline.com/ArcGIS/rest/services/Canvas/World_Light_Gray_Base/MapServer/tile/',
# 'esri-transport'='http://server.arcgisonline.com/ArcGIS/rest/services/Reference/World_Transportation/MapServer/tile/',
# 'esri-topo' = 'http://services.arcgisonline.com/ArcGIS/rest/services/World_Topo_Map/MapServer/tile/'
	
	
  

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
	crs=projection(x),
  buffer=0, fact=1,
	verbose=FALSE) {


	alltiles = osmTiles()
	pathOrig = path
	pathisname = gsub("-", ".", pathOrig) %in% 
			gsub("-", ".", names(alltiles))
	path[pathisname] = alltiles[path[pathisname]]
	
	if(length(grep("/$", path, invert=TRUE)))
		path[ grep("/$", path, invert=TRUE)] =
				paste(path[ grep("/$", path, invert=TRUE)], "/", sep="")

	if(length(grep("^http[s]*://", path, invert=TRUE)))
		path[ grep("^http[s]*://", path, invert=TRUE)] = 
				paste("http://", 
						path[ grep("^http[s]*://", path, invert=TRUE)], sep="")
	names(path) = pathOrig

	crsOut=crs
	crsIn = crs(x)
	if(all(is.na(crsIn))) {
		if(is.vector(x)){
			crsIn=crsLL
		} else{
			crsIn = crs	
		}
	}
	
	
	extMerc = .getExtent(x,crsIn, buffer, crsMerc)
  extMerc = .cropExtent(extMerc, extentMerc)
  
	if(missing(zoom)) {
	zoom = 1
	while(nTilesMerc(extMerc, zoom) <= maxTiles & zoom <= 18) {
		zoom = zoom + 1
	}
	zoom = min(c(18,max(c(1, zoom-1))))
	}
	if(verbose) cat("zoom is ", zoom, ", ", nTilesMerc(extMerc, zoom), "tiles\n")

	result = NULL
  
	for(Dpath in rev(names(path))) {
		Durl = path[Dpath]
		if(verbose){
			cat(Dpath, '\n')
			cat(Durl, '\n')
		}
		
		if(length(grep(
						'nrcan\\.gc\\.ca', Durl))){
			suffix = ''
			tileNames = 'zyx'
		} else if(length(grep('heidelberg.de/tiles/(hybrid|adminb|roadsg|roads)/?$', Durl))){
			tileNames = 'xyz='
			suffix = ''
		} else {
			suffix = '.png'
			tileNames = 'zxy'
		}
		
		thistile = try(
				getTilesMerc(extMerc, zoom=zoom,
				path=Durl,
				verbose=verbose,
				suffix=suffix,
				tileNames = tileNames),
		silent=!verbose)

		if(class(thistile)=="try-error"){
			message(paste(Durl, "not accessible"))
      thistile=NULL
		}	else {
			if(length(names(thistile))) {
				theprefix=strsplit(
						names(thistile), 
						"([rR]ed|[gG]reen|[bB]lue)$",fixed=FALSE
				)[[1]]
				names(thistile) = gsub(
						theprefix, paste(Dpath, "",sep=""), 
						names(thistile),fixed=TRUE)		
			}

		ctable = NULL
		if(!is.null(thistile)) {
			if(nlayers(thistile)==1)
				ctable = thistile@legend@colortable
			if(is.null(result)) {
				result = stack(thistile)
			} else {
	      result =  stack(result, thistile)	
			}
		}
		
	
		} # end not try-error
	} # end loop through path	

	
	if(is.null(result)) {
		# create an empty raster
		result = raster(extMerc,1,1,crs=crsMerc)
		values(result) = NA
    attributes(result)$openmap = list(
        tiles=NA,
        message=thistile,
        path=path,
        zoom=zoom
        )
	}

	
	if(!is.na(crsOut)  ){
		oldColorTable = list()
		for(D in names(result))
			oldColorTable[[D]] = result[[D]]@legend@colortable
		
	if(verbose) cat("reprojecting ", ncell(result), " cells...")
		
	# if tiles need projecting
	if(!compareCRS(projection(result), crsOut)) {

    toRaster = projectExtent(result, crsOut)
		
		# see if this raster has an unnecessarily large extent
		bboxx = try(bbox(extend(extent(x), buffer)), silent=TRUE)
		
		# check for bboxx a single point
		if(any(apply(bboxx,1,diff)<.Machine$double.eps)) {
			class(bboxx) = 'try-error'
		}
		
		projx = try(proj4string(x), silent=TRUE)
		if(class(bboxx)!="try-error" &  class(projx) != 'try-error'){
		  if(identical(projx, crsOut)) {
			  bigExtent =  extend(
						extent(x), buffer + abs(as.numeric(apply(bboxx, 1, diff))/2))
			  toRaster = raster::crop(toRaster, bigExtent)
		  }
		}
		if(any(fact > 1)){
			res(toRaster) = res(toRaster) / rep_len(fact,2)
		}

		resultProj = stack(projectRaster(result, toRaster, method="ngb"))

		for(D in names(resultProj))
			resultProj[[D]]@legend@colortable = oldColorTable[[D]]

		
    if(verbose) cat("done\n")
	} else { # crsOut and projection(result) are the same
		if(any(fact < 1)) {
			# aggregate the raster
			fact = pmax(1,round(1/fact))
			if(any(fact > 1)) {
				if(verbose) cat("aggregating tiles by ", fact,  "\n")
				resultProj = stack(aggregate(result, fact=fact, fun=min))
				for(D in names(resultProj)) resultProj[[D]]@legend = result[[D]]@legend
			}
		} else {
			if(verbose) cat("tiles arrived in the desired projection\n")
			resultProj = stack(result)
		}
	}
    
	} else { # crsOut is NA, don't project
		resultProj = stack(result)
	}


	for(D in names(resultProj)) {
    if(length(result[[D]]@legend@colortable)) {
				if(verbose) cat("copying colortable for ", D, "\n")
			resultProj[[D]]@legend@colortable =
					result[[D]]@legend@colortable
      if(any(values(resultProj[[D]])==0,na.rm=TRUE)) {
      # set NA's to transparent
      	resultProj[[D]]@legend@colortable =
          c('#FFFFFF00',
              resultProj[[D]]@legend@colortable 
          )
      	values(resultProj[[D]]) = 1+values(resultProj[[D]])
			} # end zeros
		} # end have colortable
	} # end D in names(resultProj)

	if(nlayers(resultProj)==1) 
		resultProj = resultProj[[1]]
		
  attributes(resultProj)$tiles = attributes(thistile)$tiles
  attributes(resultProj)$tiles$path = path
	
	resultProj
}

