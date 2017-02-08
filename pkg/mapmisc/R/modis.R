
getModisRaster = function() {
  modisRaster = raster(
    extent(-20015109.354,20015109.354,-10007554.677,10007554.677),
    nrow=18, ncol=36,
    crs=crsModis		
  )
  
  values(modisRaster)= 1:ncell(modisRaster)
  modisDf = data.frame(
    ID = values(modisRaster),
    h=colFromCell(modisRaster, 1:ncell(modisRaster))-1,
    v=rowFromCell(modisRaster, 1:ncell(modisRaster))-1
  )
  modisDf$hString = sprintf("%02d", modisDf$h)
  modisDf$vString = sprintf("%02d", modisDf$v)
  
  modisDf$tile = paste(
    'h', modisDf$hString, 'v', modisDf$vString, sep=''
  )
  
  modisRaster = as.factor(modisRaster)
  levels(modisRaster) = list(modisDf)
  modisRaster
}


getDegreeRaster = function() {
  degreeRaster = raster(
    extent(c(-180,180,-90,90)),
    res=1, crs=crsLL
  )
  values(degreeRaster)=1:ncell(degreeRaster)  
# coordinate is bottom left
  degreeMatXY = xyFromCell(degreeRaster, 1:ncell(degreeRaster)) - 0.5
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
  degreeRaster = ratify(degreeRaster)
  levels(degreeRaster)[[1]] = degreeDf
  degreeRaster
}

modisRaster = getModisRaster()
degreeRaster = getDegreeRaster()

getModisTiles = function(x, tiles = mapmisc::modisRaster) {
  
  xModis = projectExtent(x, projection(tiles))
  
  modisCrop = crop(tiles, 
    extend(extent(xModis), sqrt(.Machine$double.eps)), 
    snap='out')

	res = factorValues(tiles, values(modisCrop))
  
  res
}

