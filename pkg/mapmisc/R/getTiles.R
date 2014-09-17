### sample usage:
###  lancs =  getTiles(c(-2.842,-2.7579),c(54.0295,54.063),12,path="http://tile.openstreetmap.org/",maxTiles=60,verbose=TRUE)
###

#  crsMerc =CRS("+init=epsg:3857") # mercator projection
crsMerc = CRS("+proj=merc +a=6378137 +b=6378137 +lat_ts=0.0 +lon_0=0.0 +x_0=0.0 +y_0=0 +k=1.0 +units=m +nadgrids=@null +no_defs")
# crsLL = CRS("+epsg:4326")


.tile2boundingBox <- function(x,y,zoom){

  n = .tile2lat(y,zoom)
  s = .tile2lat(y+1,zoom)
  w = .tile2lon(x,zoom)
  e = .tile2lon(x+1,zoom)
  return(c(s,w,n,e))
  
}

.tile2extentMercator <- function(I, J, zoom) {
	
	# formula from http://wiki.openstreetmap.org/wiki/Slippy_map_tilenames
	n = 2 ^ zoom
	lon_rad = c(I, I+1) / n * 2*pi - pi
	lon_deg = lon_rad *180/pi
	lat_rad = atan(sinh(pi * (1 - 2 * c(J+1,J) / n)))
	lat_deg = lat_rad * 180.0 / pi

	eps=1e-8
	thePoints = raster(extent(lon_deg[1], xmax=lon_deg[2],
					ymin=lat_deg[1],ymax=lat_deg[2]), 
				crs=crsLL)
	thePointsMerc = projectExtent(thePoints, crsMerc)

	extent(thePointsMerc)			
}


.tile2lat <- function(y,zoom){
  n = pi - ((2.0 * pi * y) / (2.0^zoom))
  return ( 180.0 / pi * atan(0.5 * (exp(n) - exp(-n))))
}

.tile2lon <- function(x,zoom){
  return(((x/(2^zoom)*360)-180))
}


.lonlat2tile <- function(lon,lat,zoom){
  xtile = floor(((lon + 180) / 360) * 2^zoom)
  ytile = floor((1 - log(tan(lat*pi/180) + 1 / cos(lat*pi/180)) / pi) /2 * 2^zoom)
  return(c(xtile,ytile))
}

.getTileBounds <- function(xlim,ylim,zoom){
  LL = .lonlat2tile(xlim[1],ylim[1],zoom)
  UR = .lonlat2tile(xlim[2],ylim[2],zoom)
  return(list(LL,UR))
}

nTiles <- function(xlim,ylim,zoom){
  tb = .getTileBounds(xlim,ylim,zoom)
  nt = (tb[[2]][1]-tb[[1]][1]+1)*(tb[[1]][2]-tb[[2]][2]+1)
  if(!is.finite(nt))
	  nt = Inf
  return(nt)
}

getTilePaths <- function(xlim,ylim,zoom,path){
  tileBounds = .getTileBounds(xlim,ylim,zoom)
  LL = tileBounds[[1]]
  UR = tileBounds[[2]]
  tileData = list()
  i = 1
  for(I in LL[1]:UR[1]){
    for(J in LL[2]:UR[2]){
      tilePath = paste(path,paste(zoom,I,J,sep="/"),'.png',sep='')
      tileBounds = .tile2boundingBox(I,J,zoom)
	  tileExtent= .tile2extentMercator(I,J,zoom)
      tileData[[i]]=list(path=tilePath,bounds=tileBounds,I=I,J=J,zoom=zoom,src=path,
			  extent=tileExtent)
      i = i + 1
    }
  }
  return(tileData)
}

getTiles <- function(xlim,ylim,zoom,path,maxTiles = 16,
		cacheDir=tempdir(),
		timeOut=5*24*60*60,verbose=FALSE){
	if(verbose) {
		cat(path, "\n")
	}
	
	nt = nTiles(xlim,ylim,zoom)
  if(nt > maxTiles){
    stop("Cant get ",nt," tiles with maxTiles set to ",maxTiles)
  }
  tileData = getTilePaths(xlim,ylim,zoom,path)
  localStore = FALSE
  if(file.exists(path)){
    if(!file.info(path)$isdir){
      stop("Path ",path," is not a folder")
    }
    localStore = TRUE
  }

  thecrs = crsMerc
  
  rasters = list()
  colourtable = NULL
  for(ip in 1:length(tileData)){
   p = tileData[[ip]]$path
    if(localStore){
      where = p
    }else{
       where = .getTileCached(tileData[[ip]],cacheDir,
			   timeOut=timeOut,
			  verbose=verbose)
    }
	if(file.info(where)$size) {
		thisimage = raster::brick(where)
	} else {
		thisimage = NULL
	}
	if(!is.null(thisimage)) {

		# if only one layer
	# this must be a raster of integers referring to colour indexes	
	if(nlayers(thisimage)==1) {
		newcolours = thisimage@legend@colortable
		thisimage=thisimage[[1]]
		if(length(newcolours)) {
			if(any(is.na(values(thisimage)))) {
				newcolours[1] = NA
			}
			thisimage = thisimage + length(colourtable)
			colourtable = c(colourtable, newcolours)
		}
		names(thisimage) = gsub("^http://|/$", "", path)
	} else if (nlayers(thisimage)>1){
		thisimage = thisimage[[1:3]]

		names(thisimage) = paste(
				gsub("^http://|/$", "", path),
				c("Red","Green","Blue"),
				sep="")
	} 
	
	extent(thisimage) = tileData[[ip]]$extent
	proj4string(thisimage) = thecrs
	if(length(colourtable)) {
		thisimage@legend@colortable = colourtable
	}
	rasters[[ip]] = thisimage
	} # end not rtry error

	
	}

	if(length(rasters) > 1) {
		thenames = names(rasters[[1]])
		rasters = do.call(raster::merge, rasters)
		names(rasters) = thenames			
	} else if(length(rasters)==1){
		rasters = rasters[[1]]
	} else { # no tiles found
		rasters = NULL
	}
	if( !is.null(colourtable)) {
 		# re-order colours so most common colours are first
		newtable = unique(colourtable)
		tomatch = match(colourtable,newtable)
		rasters@legend@colortable = newtable
		values(rasters) = tomatch[values(rasters)+1]-1
	}
	
	return(rasters)	
	
}
	

.getTileCached <- function(tileData,cacheDir,timeOut=30,verbose=FALSE){
  srcDir = paste("X",make.names(tileData$src),sep="")
  tileDir = paste(file.path(cacheDir,srcDir,tileData$zoom,tileData$I))
  tilePath = paste(file.path(tileDir,tileData$J),".png",sep="")
  retrieve = FALSE
  if(!file.exists(tilePath)){
    if(verbose)cat("Getting new tile\n")
    retrieve = TRUE
  }else{
    age = difftime(Sys.time(),file.info(tilePath)$mtime,units = "mins")
    if (age > timeOut){
      if(verbose)cat("Tile aged ",age," expired from cache\n")
      retrieve = TRUE
    }else if(file.info(tilePath)$size<1) {
		if(verbose)cat("Tile file is too small, retreiving again\n")
		retrieve = TRUE
	} else {
      if(verbose){cat("Tile found in cache\n")}
    }
  }
  if(retrieve){
    if(verbose)cat("Downloading tile\n")
    dir.create(tileDir,recursive=TRUE,showWarnings=FALSE)
    p = tileData$path
    res=try(download.file(p,tilePath,mode="wb",quiet=!verbose),silent=TRUE)
	if(class(res)=="try=error"){
		message(p, "not accessible")
	}
}
  
  return(tilePath)
}
