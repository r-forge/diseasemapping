
tonerToTrans = function(x, pattern= "(red|green|blue)$", 
  power=0.5, col='black',
  threshold=Inf) {

  maxColorValue = max(maxValue(x))
  tilesOrig = attributes(x)$tiles
  
  rgbLayers = grep(pattern, 
    names(x), 
    ignore.case=TRUE)
  
  if(!length(rgbLayers))
    warning("can't find specified pattern in x")
  
  if(length(rgbLayers)> 1) {
    xMax = calc(x[[rgbLayers]], 
      min, na.rm=FALSE)
  } else {
    xMax = x[[rgbLayers]]
  }
  if(threshold < maxColorValue)
    xMax = reclassify(xMax, 
      cbind(threshold, maxColorValue, maxColorValue))
  
  newTrans = ( 
      (maxColorValue-values(xMax) 
        ) / maxColorValue)^power
  newTrans[is.na(newTrans)] = 0
  
  col = col2rgb(col)[,1]/255
  
  newCol = grDevices::rgb(
    col['red'],
    col['green'],
    col['blue'],
    newTrans
  )
  newCol = factor(newCol)
  levels(newCol) = gsub("FFFFFFFF", "FFFFFF00", levels(newCol))
  result = raster(x)
  values(result) = as.numeric(newCol)
  result@legend@colortable = c(NA,levels(newCol))

  names(result) = gsub("(red|green|blue)$", "", 
    names(x)[rgbLayers[1]], ignore.case=TRUE)
  
  attributes(result)$tiles = tilesOrig
  attributes(result)$tiles$tonerToTrans = match.call()
  result
}


rgbtToIndex = function(x, pattern="(red|green|blue|trans)$") {
  
  rgbLayers = grep(pattern, 
    names(x), 
    ignore.case=TRUE)
  
  
  if(length(rgbLayers) != 4)
    warning("x doesn't seem to be RGB trans")
  
  xsub = subs(
    x[[rgbLayers]], 
    data.frame(NA,0), 
    subsWithNA=FALSE)
  
  
  xMaxValue = max(maxValue(x), na.rm=TRUE)
  
  newCol = grDevices::rgb(
    values(xsub[[grep("red$", names(xsub), ignore.case=TRUE)[1]]]),
    values(xsub[[grep("green$", names(xsub), ignore.case=TRUE)[1]]]),
    values(xsub[[grep("blue$", names(xsub), ignore.case=TRUE)[1]]]),
    values(xsub[[grep("trans$", names(xsub), ignore.case=TRUE)[1]]]),
    maxColorValue = xMaxValue
  )
  newCol = factor(newCol)		
  
  newCol[grep("00$", newCol)] = NA
  
  newCol = factor(newCol)
  
  result = raster(x)
  
  names(result) = gsub(
    "(red|green|blue|trans)$", "", names(x), 
    ignore.case=TRUE)[1]
  
  values(result) = as.numeric(newCol)
  
  result@legend@colortable = c(NA,levels(newCol))
  
  result
}
