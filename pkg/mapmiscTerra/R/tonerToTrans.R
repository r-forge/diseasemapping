
tonerToTrans = function(x, pattern= "(red|green|blue)$", 
  power=0.5, col='black',
  threshold=Inf) {

  # the colorize function doesn't appear to work
  maxColorValue = max(terra::minmax(x))
  xValues = terra::values(x)
  xUnique = xValues[!duplicated(xValues),]
  xUnique = cbind(value = seq(0, nrow(xUnique)-1), xUnique)

  xMax = apply(xUnique[,2:4], 1, min)
  if(threshold < maxColorValue)
    xMax[which(xMax > threshold)] = maxColorValue 

  newTrans = floor(maxColorValue * ( 
      (maxColorValue-xMax)  / maxColorValue)^power)
  newTrans[is.na(newTrans)] = 0
  xUnique[,5] = newTrans


  xValues = cbind(cell = 1:nrow(xValues), xValues[,1:3])

  result = rast(x, nlyrs=1)

  valuesOfResult = merge(xValues, xUnique, all.x=TRUE, by.x=2:4, by.y=2:4)
  valuesOfResult = valuesOfResult[match(1:nrow(valuesOfResult), valuesOfResult[,'cell']), ]

  terra::values(result) = valuesOfResult[,'value']
  terra::coltab(result) = xUnique

  attributes(result)$tiles = attributes(x)$tiles
  attributes(result)$openmap = attributes(x)$openmap
  
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
