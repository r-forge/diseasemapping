
tonerToTrans = function(x, pattern= "(red|green|blue)$", 
  power=0.5, col='black',
  threshold=Inf) {

  # the colorize function doesn't appear to work
  maxColorValue = max(terra::minmax(x))
  xValues = terra::values(x)
  if(!all(c('red','green','blue') %in% colnames(xValues))) {
    warning("x needs red, green, blue layers")
  }
  xUnique = xValues[!duplicated(xValues[,c('red','green','blue')]),]
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

  valuesOfResult = merge(as.data.frame(xValues), as.data.frame(xUnique), by = c('red','green','blue'))
  valuesOfResult = valuesOfResult[match(1:nrow(valuesOfResult), valuesOfResult[,'cell']), ]

  terra::values(result) = valuesOfResult[,'value']
  terra::coltab(result) = xUnique

  attributes(result)$tiles = attributes(x)$tiles
  attributes(result)$openmap = attributes(x)$openmap
  
  attributes(result)$tiles$tonerToTrans = match.call()
  result
}

