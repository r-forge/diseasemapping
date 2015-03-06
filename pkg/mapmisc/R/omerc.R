omercProj4string = function(
    lon, lat, angle, 
    x=0,y=0, inverseAngle=0,
    scale=1,
    ellps='WGS84', units='m') {
  
  negAngle = angle<0
  angle[negAngle] = 360 - angle[negAngle]
  angle[angle==90]=89

  whichZeros = angle==0
  
  result = paste(
      "+proj=omerc",
      " +lat_0=", lat,
      " +lonc=", lon,
      " +gamma=", inverseAngle,
      " +k=", scale, 
      " +alpha=", angle, 
      " +x_0=", x,
      " +y_0=", y,
      " +ellps=", ellps,
      " +units=", units,
      sep="")
  
  if(any(whichZeros)) {
    result[whichZeros] = 
        gsub("omerc","tmerc", result[whichZeros])
    result[whichZeros] = 
        gsub(" \\+(alpha|gamma)=([[:digit:]]|\\.)+","", 
            result[whichZeros])
    result[whichZeros] = 
        gsub("lonc=","lon_0=", 
            result[whichZeros])
    
  }
  
  result
}

omerc = function(
    x, angle=0, undo=FALSE
) {
  
  
  crs = projection(x)
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

  if(undo){
    inverseAngle = -angle
  } else {
    inverseAngle = 0
  }
  
  rotatedProj4string = 
      omercProj4string(
      lon=theCentre[1],
      lat=theCentre[2], 
      angle=angle,
      inverseAngle = inverseAngle
      )

  rotatedCRS = lapply(rotatedProj4string, CRS)
  
  if(is.numeric(x) | !haveRgdal) {
      if(length(rotatedCRS)==1)
        rotatedCRS = rotatedCRS[[1]]
      return(rotatedCRS)
  }

  # if we're undoing the rotation
  # find the optimal rotatino
  # for a small bounding box
  # before the rotation is undone
  if(undo){
    rp2 = gsub("gamma=([[:digit:]]|\\.)+", 
        "gamma=0", rotatedProj4string)
    rc2 = lapply(rp2, CRS)
  } else {
    rc2 = rotatedCRS
  }
    
    xTrans = mapply(
        function(CRSobj) {
          abs(prod(apply(bbox(
               spTransform(x, CRSobj)           
                          ), 1, diff)
              ))
        },
        CRSobj=rc2
        )
    rotatedCRS = rotatedCRS[[
        which.min(xTrans)
        ]]
        
  result = spTransform(x, rotatedCRS)
        
  result
  
}
  
  