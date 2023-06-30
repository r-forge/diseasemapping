
getRowCol <- function(extMerc,
  zoom, 
  rasterSphere = .getRasterMerc(zoom)){
  
  Sy=terra::rowFromY(rasterSphere, c(
      terra::ymax(extMerc), terra::ymin(extMerc)
    ))
  Sx=terra::colFromX(rasterSphere, c(
      terra::xmin(extMerc), terra::xmax(extMerc)
    ))
  
  
  list(
    col = seq(Sx[1],Sx[2]),
    row = seq(Sy[1],Sy[2])
  ) 
}


getRasterNrcanDontNeed = function(zoom) {
  
  # raster for maps from 
  # http://geoappext.nrcan.gc.ca/arcgis/rest/services/BaseMaps/CBMT_CBCT_GEOM_3978/MapServer
  
  origin = c(-3.46558E7, 3.931E7) 
  nrCrs = terra::crs(
  "+proj=lcc +lat_1=49 +lat_2=77 +lat_0=49 +lon_0=-95 +x_0=0 +y_0=0 +ellps=GRS80 +towgs84=0,0,0,0,0,0,0 +units=m +no_defs") 
  
  
  nrRes = c(38364.660062653464,
    22489.62831258996,
    13229.193125052918,
    7937.5158750317505,
    4630.2175937685215,
    2645.8386250105837,
    1587.5031750063501,
    926.0435187537042,
    529.1677250021168,
    317.50063500127004,
    185.20870375074085,
    111.12522225044451,
    66.1459656252646,
    38.36466006265346,
    22.48962831258996,
    13.229193125052918,
    7.9375158750317505,
    4.6302175937685215
  )
  
  eps=0.1
  worldLL = terra::rast(
    terra::ext(-180+eps,180-eps,-90+eps,90-eps), 
    crs=crsLL,
    res=0.01
  )
  
  worldNrcan = terra::project(worldLL, nrCrs, res=1000)
  
  nrcanExtent = terra::ext(
    -4282638.06150141,
    4852210.1755664,
    -5153821.09213678,
    4659267.000000001)
  
  nrcanExtent = terra::ext(
    -7786476.885838887,
    7148753.233541353,
    -927807.6591668031,
    7928343.534071138
  )
  
  nrcanRaster = terra::rast(
    nrcanExtent,
    nrow = 2^(1+zoom),
    ncol = 2^(1+zoom),
    crs = nrCrs
  )		
  
  N = 2^zoom 
  
  terra::rast(extentMerc, nrows = N, ncols=N, crs=crsMerc)
}

getRowColNrcanDontNeed <- function(
  extMerc,
  zoom, 
  rasterSphere = .getRasterMerc(zoom)){
  
  
  Sy=terra::rowFromY(rasterSphere, c(
      terra::ymax(extMerc), terra::ymin(extMerc)
    ))
  Sx=terra::colFromX(rasterSphere, c(
      terra::xmin(extMerc), terra::xmax(extMerc)
    ))
  
  
  list(
    col = seq(Sx[1],Sx[2]),
    row = seq(Sy[1],Sy[2])
  ) 
}


nTilesMerc <- function(extMerc,zoom){
  
  SrowCol = terra::getRowCol(extMerc, zoom=zoom)
  
  length(SrowCol[[1]])*length(SrowCol[[2]])
  
}


getTilesMerc = function(
  extMerc=extentMerc, 
  zoom=1, 
  path="http://tile.openstreetmap.org/",
  cachePath='.',
  cacheDir= make.names(gsub(
          "^http.*//([[:alpha:]][.])*((tile|basemap)s?[.][[:digit:]]?)?(openstreetmap[.])?|[[:punct:]]$", 
          "", path)),
  verbose=FALSE, suffix = '.png',
  tileNames = 'zxy'){
  
  
  cacheDir = file.path(cachePath, cacheDir)
  
  rasterSphere = .getRasterMerc(zoom)  
  
  SrowCol = getRowCol(extMerc, rasterSphere=rasterSphere)
  
  rasters = list()
  colourtable = NULL
  
  
  for(Dx in SrowCol$col){
    
    rasterCol= list()
    for(Dy in SrowCol$row) {
      
      Dcell = terra::cellFromRowCol(rasterSphere, Dy, Dx)

      # extent of the cell
      Dextent = terra::ext(rep(terra::xyFromCell(rasterSphere, Dcell), each=2) + rep(terra::res(rasterSphere), each=2)*c(-1,1,-1,1)/2)
      
      
      if(tileNames == 'zxy') {
        Dcache = file.path(cacheDir, zoom, Dx-1)
        Dpath = paste(path,zoom,'/',Dx-1,'/',sep='')
        Dtile = paste(Dy-1, suffix,sep='')
        Durl = paste(Dpath, Dtile, sep='')
      } else if (tileNames == 'zxy=') {
        Dcache = cacheDir
        Dtile = paste('z=',zoom, '&x=', Dx-1, '&y=',Dy-1,suffix, sep='')
        Durl = paste(path, Dtile, sep='')				
      } else if (tileNames == 'xyz=') {
        Dcache = cacheDir
        Dtile = paste('x=',Dx-1, '&y=', Dy-1, '&z=',zoom,suffix, sep='')
        Durl = paste(path, Dtile, sep='')				
      } else if (tileNames == 'xyz') {
        Dcache = cacheDir
        Dtile = paste(Dx-1, '/', Dy-1, '/',zoom,suffix, sep='')
        Durl = paste(path, Dtile, sep='')				
      } else if (tileNames == 'zyx') {
        Dcache = file.path(cacheDir, zoom, Dy-1)
        Dpath = paste(path,zoom,'/',Dy-1,'/',sep='')
        Dtile = paste(Dx-1, suffix, sep='')
        Durl = paste(Dpath, Dtile, sep='')
      } else {
        warning('tileNames must be zxy or zyx or xyz=')
      }
      
      dirCreateMapmisc(path=Dcache,recursive=TRUE,showWarnings=FALSE)
      Dfile = file.path(Dcache, Dtile)
      
      
      Dsize = file.info(Dfile)['size']
      if(!any(Dsize > 0,na.rm=TRUE)) {
        if(verbose) cat("downloading ", Durl, "\n")
        try(downloadFileMapmisc(
            Durl, Dfile, quiet=!verbose, 
            method='auto', mode = 'wb'
          ), silent=TRUE)
      } else {
        if(verbose) cat("tile ", Dfile, " cached\n")
      }
      
      thisimage = try(terra::rast(Dfile), silent=TRUE)
      
      if(any(class(thisimage)=='try-error')) {
        if(verbose) warning("tile ", Dfile, " cannot be loaded")
        
        thisimage = terra::rast(
          Dextent, nrows=256, ncols=256, crs=crsMerc
        )
        terra::values(thisimage) = NA
      } else {
        terra::crs(thisimage) = crsMerc
        terra::ext(thisimage) = Dextent
      }
      
      if(terra::nlyr(thisimage)==1) {
        # single layer, there's a colortable
        newcolours = terra::coltab(thisimage)
        if(length(newcolours)) {
          if(any(is.na(values(thisimage)))) {
            newcolours[1] = NA
          }
          thisimage = thisimage + length(colourtable)
          colourtable = c(colourtable, newcolours)
          terra::coltab(thieimage) = colourtable
        }
        names(thisimage) = gsub("^http://|/$", "", path)
      } else if (terra::nlyr(thisimage)>1){
        
        # if values are 1 to 256, change to 0 to 255
        terra::setMinMax(thisimage)
        if(max(terra::minmax(thisimage))==256 ) {
          thisimage = thisimage - 1
        }
        
        cnames = c('Red','Green','Blue','Trans')[1:min(c(4,terra::nlyr(thisimage)))]
        
        names(thisimage)[1:length(cnames)] = paste(
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
      } # end more than one layer 
      
      rasterCol[[paste('y', Dy,sep='')]] = thisimage
    } # end Dy
    
    if(length(rasterCol)> 1) {
      names(rasterCol)[1:2] = c('x','y')
      rasters[[paste('x',Dx,sep='')]] = do.call(
        raster::merge, rasterCol
      )
      names(rasters[[paste('x',Dx,sep='')]]) = 
        names( rasterCol[[1]])
    } else {
      rasters[[paste('x',Dx,sep='')]] = rasterCol[[1]]
    }
  } # end Dx
  
  # merge the rows
  if(length(rasters) > 1) {
    thenames = names(rasters[[1]])
    
    names(rasters)[1:2] = c('x','y')
    rastersMerged = do.call(
      raster::merge, rasters
    )
    names(rastersMerged) = names(rasters[[1]])
  } else {
    rastersMerged = rasters[[1]]
  }
  
  # add the colortable
  if( length(colourtable)) {
    # re-order colours so most common colours are first
    newtable = unique(colourtable)
    tomatch = match(colourtable,newtable)
    rastersMerged@legend@colortable = newtable
    values(rastersMerged) = tomatch[values(rastersMerged)+1]-1
  }
  attributes(rastersMerged)$tiles = 
    list(tiles = length(SrowCol[[1]])*length(SrowCol[[2]]), 
      zoom=zoom,
      path=path)
  
  return(rastersMerged)	
  
}


