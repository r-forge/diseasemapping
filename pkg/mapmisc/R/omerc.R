
omercCrs = function(lat, lon, angle) {
  
  if(angle<0)
    angle = 360-angle
  if(angle==90)
    angle = 89
  if(angle==0){
    result = CRS(paste(
            "+proj=tmerc +lat_0=",
        lat,
        " +lon_0=",
        lon,
            " +k=1",
            " +x_0=0 +y_0=0 +ellps=WGS84 +units=m",
            sep=""))
  } else {
    result = CRS(paste(
        "+proj=omerc +lat_0=",
        lat,
        " +lonc=",
        lon,
        " +gamma=0.0 +k=1 +alpha=",
        angle, 
        " +x_0=0 +y_0=0 +ellps=WGS84 +units=m",
          sep=""))
  }	 
  result
}

omerc = function(
    x, angle=0, crs=projection(x)
) {
  
  
  if(is.na(crs)){
    crs = crsLL
  }
  
  if(is.character(crs))
    crs = CRS(crs)
  
  haveRgdal = requireNamespace('rgdal', quietly=TRUE)

  if(!is.numeric(x)){
    theCentre = bbox(x)
    theCentre = theCentre[,'min'] +
      apply(theCentre, 1, diff)/2
  } else {
    theCentre = x[1:2]
  }
  
  if(!isLonLat(crs)) {
      theCentre = SpatialPoints(
          t(theCentre[1:2]),
          proj4string=crs
          )
       if(!haveRgdal) {
          warning('install rgdal to use projections other than long-lat')
       } else {
        theCentre = spTransform(theCentre, crsLL)
      }
      theCentre = as.vector(theCentre@coords)
  } # crs not LL

  rotatedCrs = mapply(
      omercCrs,
      angle=angle,
      MoreArgs=list(
      lat=theCentre[2], lon=theCentre[1])
  )
  
    if(is.numeric(x) | !haveRgdal) {
      if(length(rotatedCrs)==1)
        rotatedCrs = rotatedCrs[[1]]
      return(rotatedCrs)
    }
    
    xTrans = mapply(
        spTransform,
        CRSobj=rotatedCrs,
        MoreArgs=list(x=x)
        )
    
    bbArea = mapply(
        function(stuff){
          abs(prod(apply(bbox(stuff),1,diff)))
        },
        stuff=xTrans
        )
    
    return(xTrans[[which.min(bbArea)]])
    
}
  
  