tpeqdProj4string = function(
    lon1, lat1, lon2, lat2, 
    x=0,y=0, datum='WGS84',
    ellps='WGS84', units='m',
    axis='enu',
    crs=TRUE) {
  

  res = paste(
      "+proj=tpeqd",
      " +lat_1=",lat1,
      " +lon_1=",lon1,
      " +lat_2=",lat2,
      " +lon_2=",lon2,
      " +x_0=",x,
      " +y_0=",y,
      " +axis=", axis,
      " +ellps=",ellps,
#      " +datum=",datum,
      " +units=",units,
      " +no_defs",sep=''
      )
  if(crs){
    res = lapply(res, crs)
  }
  res
}

tpeqd = function(x, offset=c(0,0), axis='enu', ...){
  
  if(length(grep("^SpatVector", class(x)))){
    x = project(x, crsLL)
    x = terra::crds(x)
  }
  
  x = as.matrix(x)[1:2,1:2]
#  x = x[order(x[,2],decreasing=TRUE),]

    # check if crossing -180
    if(x[1,1] > x[2,1]) {
      x[1,1] = x[1,1] - 360
    }

  
  
  res= tpeqdProj4string(
      x[1,1],x[1,2],x[2,1],x[2,2],
      x=offset[1],y=offset[2],
      axis=axis, crs=FALSE
      )
  if(length(res)[[1]])
    res = res[[1]]
      
  circle = exp(1i*seq(0, 2*pi, len=1001)[-1])
  circle = cbind(Re(circle), Im(circle))
  offsetmat = matrix(offset, nrow(circle), 2, byrow=TRUE)
  step = 100*1000
  radiusScale = ceiling(extentMerc[2]/step)*step
  somePointsOut = TRUE
  while( (radiusScale > 0) & somePointsOut) {
    radiusScale = radiusScale - step
    xxOrig = vect(radiusScale * circle + offsetmat, crs=res, type='points')
    xxLL = suppressWarnings(project(xxOrig, crsLL))
    somePointsOut = any(is.nan(as.vector(geom(xxLL)[,c('x', 'y')])))
  }
  step = round(step/10)
  # try to make one axis bigger 
  radiusScale = rep(radiusScale, 2)
  noPointsOut = TRUE
  while(noPointsOut & all(radiusScale < extentMerc[2])) {
    radiusScale = radiusScale + c(step,0)
    xxOrig = vect(matrix(radiusScale, nrow(circle), 2, byrow=TRUE) * circle + offsetmat, 
      crs=res, type='points')
    xxLL = suppressWarnings(project(xxOrig, crsLL))
    noPointsOut = !any(is.nan(as.vector(geom(xxLL)[,c('x', 'y')])))
  }
  radiusScale = radiusScale - c(step,0)
  noPointsOut = TRUE
  while(noPointsOut & all(radiusScale < extentMerc[2])) {
    radiusScale = radiusScale + c(0,step)
    xxOrig = vect(matrix(radiusScale, nrow(circle), 2, byrow=TRUE) * circle + offsetmat, 
      crs=res, type='points')
    xxLL = suppressWarnings(project(xxOrig, crsLL))
    noPointsOut = !any(is.nan(as.vector(geom(xxLL)[,c('x', 'y')])))
  }
  radiusScale = radiusScale - c(0,step)


  theEllipse = vect(matrix(radiusScale, nrow(circle), 2, byrow=TRUE) * circle + offsetmat, crs=res, type='polygons')
  thebox = suppressWarnings(llCropBox(res, ellipse = theEllipse, ...))
  res = terra::crs(res)
   
  attributes(res)[names(thebox)] = thebox

  res
  
}