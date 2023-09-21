
tonerToTrans = function(x, pattern= "(red|green|blue)$", 
  power=0.5, col='black',
  threshold=Inf) {

  if(all(terra::has.colors(x))) {
    xValues = terra::coltab(x)[[1]]    
    result = terra::deepcopy(x)
  } else {
    xValues = cbind(value=1:terra::ncell(x), terra::values(x))
    result = rast(x, nlyrs=1)
    terra::values(result) = 1:terra::ncell(result)
  }
  Scol = grep(pattern, colnames(xValues), value=TRUE)
  if(!length(Scol)) {
    warning("cant find RGB columns")
    Scol = grep("value$", colnames(xValues), invert=TRUE, value=TRUE)
  }
  maxColorValue = max(xValues[,Scol])

  xMax = apply(xValues[,Scol,drop=FALSE], 1, min)
  if(threshold < maxColorValue)
    xMax[which(xMax > threshold)] = maxColorValue 

  newTrans = floor(maxColorValue * ( 
      (maxColorValue-xMax)  / maxColorValue)^power)
  newTrans[is.na(newTrans)] = 0
  if(!any(colnames(xValues)=='alpha')) xValues = cbind(xValues, alpha = NA)
  xValues[,'alpha'] = newTrans

  xUnique = unique(newTrans)
  xUnique = cbind(value = 1:length(xUnique), alpha = xUnique)

  if(all(terra::has.colors(x))) {
    xValues = cbind(xValues,
     newvalue = xUnique[match(xValues[,'alpha'], xUnique[,'alpha']), 'value'])
    terra::values(result) = xValues[match(terra::values(result), xValues[,'value']), 'newvalue']
  } else {
    terra::values(result) = xUnique[match(xValues[,'alpha'], xUnique[,'alpha']), 'value']
  }

  if(is.character(col)) {
      col = drop(grDevices::col2rgb(col[1]))
  }
  theColtab = cbind(
      value = xUnique[,'value'],
      matrix(col, nrow(xUnique),  length(col), byrow=TRUE,
        dimnames = list(NULL, names(col))),
      alpha=xUnique[,'alpha']
    )

  terra::coltab(result) = theColtab

  attributes(result)$tiles = attributes(x)$tiles
  attributes(result)$openmap = attributes(x)$openmap
  
  attributes(result)$tiles$tonerToTrans = match.call()
  result
}

