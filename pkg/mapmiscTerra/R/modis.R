
getModisRaster = function() {
  modisRaster = rast(
    terra::ext(-20015109.354,20015109.354,-10007554.677,10007554.677),
    nrow=18, ncol=36,
    crs=crsModis		
  )
  
  terra::values(modisRaster)= 1:terra::ncell(modisRaster)
  modisDf = data.frame(
    ID = terra::values(modisRaster),
    h=terra::colFromCell(modisRaster, 1:terra::ncell(modisRaster))-1,
    v=terra::rowFromCell(modisRaster, 1:terra::ncell(modisRaster))-1
  )
  modisDf$hString = sprintf("%02d", modisDf$h)
  modisDf$vString = sprintf("%02d", modisDf$v)
  
  modisDf$tile = paste(
    'h', modisDf$hString, 'v', modisDf$vString, sep=''
  )
  
  modisRaster = terra::as.factor(modisRaster)
  levels(modisRaster) = list(modisDf)
  modisRaster
}


getDegreeRaster = function() {
  degreeRaster = rast(
    terra::ext(c(-180,180,-90,90)),
    res=1, crs=crsLL
  )
  terra::values(degreeRaster)=1:terra::ncell(degreeRaster)  
# coordinate is bottom left
  degreeMatXY = terra::xyFromCell(degreeRaster, 1:terra::ncell(degreeRaster)) - 0.5
  degreeDf = data.frame(
    ID=1:nrow(degreeMatXY),
    x=degreeMatXY[,'x'],
    y=degreeMatXY[,'y'],
    ns = paste(c('n','s')[1+(degreeMatXY[,'y'] < 0 ) ], 
      sprintf("%03d", abs(degreeMatXY[,'y']) ), sep=''
    ),
    ew = paste(c('e','w')[1+(degreeMatXY[,'x'] < 0 ) ], 
      sprintf("%03d", abs(degreeMatXY[,'x']) ), sep=''
    ),
    stringsAsFactors=FALSE
  )
  degreeDf$tile = paste(degreeDf$ns, degreeDf$ew, sep='')
  levels(degreeRaster) = degreeDf
  degreeRaster
}


getModisTiles = function(x, tiles) {
 
  if(missing(tiles)) tiles = getDegreeRaster()

  xModis = project(x, crs(tiles))
  
  modisCrop = terra::crop(tiles, 
    terra::extend(terra::ext(xModis), sqrt(.Machine$double.eps)), 
    snap='out')

	res = terra::cats(tiles)[[1]]
  names(res) = gsub("lyr.1", "ID", names(res))
  res = res[match(terra::values(modisCrop), res$ID), ]
  
  res
}

