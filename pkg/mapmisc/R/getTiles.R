
openmapExtentLL = extent(-180, 180,-85.05113,85.05113)
crsSphere = CRS("+proj=longlat +ellps=sphere")
crsMercSphere = CRS("+proj=merc +ellps=sphere +units=m")

if(FALSE) {
  worldLimsLL = SpatialPoints(
      t(bbox(openmapExtentLL)), 
      proj4string=crsLL)
  as.vector(extent(spTransform(worldLimsLL, crsMerc)))
  as.vector(extent(spTransform(worldLimsLL, crsMercSphere)))
} 
openmapExtentMercSphere = extent(-20015077,  20015077, -20015079,  20015079)

getRasterSphere = function(zoom) {
  N = 2^zoom 
  raster(openmapExtentMercSphere, nrow = N, ncol=N, crs=crsMercSphere)
}

.tile2extentMercator <- function(I, J, zoom) {
  
  rastSphere = getRasterSphere(zoom)
  
  extent(
      c(xFromCol(rastSphere, I+1) + c(-1,1)*xres(rastSphere)/2,
          yFromRow(rastSphere, J+1) + c(-1,1)*yres(rastSphere)/2))
}

tileEps = sqrt(.Machine$double.neg.eps)

.getTileBounds <- function(xlim,ylim,zoom){
  LL = .lonlat2tile(xlim[1],ylim[1],zoom)
  UR = .lonlat2tile(
      lon=xlim[2]- tileEps,
      lat=ylim[2]- tileEps,
      zoom)
  return(list(LL,UR))
}

.lonlat2tile <- function(lon,lat,zoom){
  xtile = pmax(0,floor(((lon + 180) / 360) * 2^zoom))
  ytile = pmax(0,floor((1 - log(tan(lat*pi/180) + 1 / cos(lat*pi/180)) / pi) /2 * 2^zoom))
  return(c(xtile,ytile))
}


getTileRowCol <- function(extLL,zoom){
  
  
  xlim= c(xmin(extLL), xmax(extLL))
  ylim = c(ymin(extLL), ymax(extLL))
  
  
  forS = .getTileBounds(xlim,ylim,zoom)
  
  Sxy = list(
      col = sort(seq(forS[[1]][1], forS[[2]][1])),
      row = sort(seq(forS[[1]][2], forS[[2]][2]))
  )
  
  return(Sxy)
}




openmapExtentMerc = extent(-20037508,  20037508, -19994877,  19994877)
getRasterMerc = function(zoom) {
  N = 2^zoom 
  raster(openmapExtentMerc, nrow = N, ncol=N, crs=crsMerc)
}

nTiles <- function(extLL,zoom){
  
  prod(unlist(lapply(getTileRowCol(extLL, zoom), length)))
  
}


.tile2boundingBox <- function(x,y,zoom){

  n = .tile2lat(y,zoom)
  s = .tile2lat(y+1,zoom)
  w = .tile2lon(x,zoom)
  e = .tile2lon(x+1,zoom)
  return(c(s,w,n,e))
  
}


.tile2lat <- function(y,zoom){
  n = pi - ((2.0 * pi * y) / (2.0^zoom))
  return ( 180.0 / pi * atan(0.5 * (exp(n) - exp(-n))))
}

.tile2lon <- function(x,zoom){
  return(((x/(2^zoom)*360)-180))
}



 getTiles <- function(extLL,zoom,path,
		cacheDir=tempdir(),
		timeOut=5*24*60*60,verbose=FALSE){
	if(verbose) {
		cat(path, "\n")
	}
  
  tileData = getTileRowCol(extLL,zoom)

  localStore = FALSE
  if(file.exists(path)){
    if(!file.info(path)$isdir){
      stop("Path ",path," is not a folder")
    }
    localStore = TRUE
  }

  
  rasters = list()
  colourtable = NULL
    
  for(Dy in tileData$row) { 
    rasterRow = list()
    for(Dx in tileData$col){
       
     p = paste(path,paste(zoom,Dx,Dy,sep="/"),'.png',sep='')
   
      if(localStore){
        where = p
      } else {
      
       where = .getTileCached(
           src=path, 
           path=p,
           col=Dx, row=Dy,
           zoom=zoom,
           cacheDir=cacheDir,
			   timeOut=timeOut,
			  verbose=verbose)
      }
	
      thisimage = try( raster::brick(where), silent=TRUE)

      if(class(thisimage)=='try-error') {
        if(verbose) warning("tile ", path, " cannot be loaded")
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

    cnames = c('Red','Green','Blue','Trans')[1:min(c(4,nlayers(thisimage)))]
    
		names(thisimage) = paste(
				gsub("^http://|/$", "", path),
				cnames,
				sep="")
    
    transLayer = grep("Trans$", names(thisimage), value=TRUE)
    colLayer = grep("Trans$", names(thisimage), value=TRUE,invert=TRUE)
    if(length(transLayer)) {
      # convert transparent to NA
      transLayer = reclassify(thisimage[[transLayer]], data.frame(-Inf, 10, NA))  
      thisimage = raster::mask(thisimage, transLayer)
    }
	} # end no colortable 
	
	extent(thisimage) = .tile2extentMercator(Dx, Dy, zoom) 
	proj4string(thisimage) = crsMercSphere
	if(length(colourtable)) {
		thisimage@legend@colortable = colourtable
	}
	rasterRow[[paste('x', Dx,sep='')]] = thisimage
	} # end not rtry error
  } # end Dx
  if(length(rasterRow)> 1) {
    names(rasterRow)[1:2] = c('x','y')
	  rasters[[paste('y',Dy,sep='')]] = do.call(
      raster::merge, rasterRow
          )
   } else {
     rasters[[paste('y',Dy,sep='')]] = rasterRow[[1]]
   }
	} # end Dy

  # merge the rows
  if(length(rasters) > 1) {
    thenames = names(rasters[[1]])
    
    names(rasters)[1:2] = c('x','y')
    rastersMerged = try(do.call(
      raster::merge, rasters
    ), silent=TRUE)

  if(class(rastersMerged)=='try-error'){
  # probably near the poles.  rasters have different y res
  
  baseRow = which.max(unlist(lapply(rasters, yres)))
  newExtent = extent(rasters[[baseRow]])
  ymax(newExtent) = max(unlist(lapply(rasters, ymax)))
  ymin(newExtent) = min(unlist(lapply(rasters, ymin)))
  
  resTemplate = extend(raster(rasters[[baseRow]]), newExtent)
  
  for(Dy in names(rasters)) {
  rasters[[Dy]] = projectRaster(rasters[[Dy]],  
      to=crop(resTemplate, extent(rasters[[Dy]])),
      method='ngb'
  )
  } 
  rasters = do.call(raster::merge, rasters)
  names(rasters) = thenames
  } else { # end try error
    rasters = rastersMerged
  }

  } else { # end length(rasters)>1
    rasters = rasters[[1]]
  }
  
  # add the colortable
  if( !is.null(colourtable)) {
 		# re-order colours so most common colours are first
		newtable = unique(colourtable)
		tomatch = match(colourtable,newtable)
		rasters@legend@colortable = newtable
		values(rasters) = tomatch[values(rasters)+1]-1
	}
  attributes(rasters)$tiles = 
      list(tiles = length(tileData), 
          zoom=zoom,
          path=path)
	
	return(rasters)	
	
}
	

.getTileCached <- function(src, path, row, col, zoom, cacheDir,timeOut=30,verbose=FALSE){
  srcDir = paste("X",make.names(src),sep="")
  tileDir = paste(file.path(cacheDir,srcDir,zoom,col))
  tilePath = paste(file.path(tileDir,row),".png",sep="")
  retrieve = FALSE
  if(!all(file.exists(tilePath))){
    if(verbose)cat("Getting new tile\n")
    retrieve = TRUE
  }else{
    age = difftime(Sys.time(),file.info(tilePath)$mtime,units = "mins")
    if (any(age > timeOut)){
      if(verbose)cat("Tile aged ",age," expired from cache\n")
      retrieve = TRUE
    } else if(any(file.info(tilePath)$size<1)) {
		  if(verbose)cat("Tile file is too small, retreiving again\n")
		    retrieve = TRUE
	  } else {
      if(verbose){cat("Tile found in cache\n")}
    }
  }
  if(retrieve){
    if(verbose)cat("Downloading tile\n")
    dir.create(tileDir,recursive=TRUE,showWarnings=FALSE)
    for(p in path) {
      res=try(download.file(p,tilePath,mode="wb",quiet=!verbose),silent=TRUE)
	    if(class(res)=="try=error"){
		    message(p, "not accessible")
  	  }
    }
  }
  
  return(tilePath)
}
