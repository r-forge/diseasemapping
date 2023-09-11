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

tpeqd = function(x, offset=c(0,0), axis='enu'){
  
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
      
  theEllipse = vect(crsRegionEllipse(res, offset), crs=res, type='polygons')

  thebox = suppressWarnings(llCropBox(res, ellipse = theEllipse, remove.holes=TRUE, crop.poles=TRUE))
  res = terra::crs(res)
   
  attributes(res)[names(thebox)] = thebox

  res
  
}