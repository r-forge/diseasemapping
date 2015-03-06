omercProj4string = function(
    lon, lat, angle, 
    x=0,y=0, inverseAngle=0,
    scale=1,
    ellps='WGS84', units='m') {
  
  negAngle = angle<0
  angle[negAngle] = 360 + angle[negAngle]
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

  
  rotatedProj4string = 
      omercProj4string(
      lon=theCentre[1],
      lat=theCentre[2], 
      angle=angle)

  rotatedCRS = lapply(rotatedProj4string, CRS)
  
  if(undo & haveRgdal) {
    pointNorth = SpatialPoints(
        rbind(
            theCentre,
            theCentre + c(0, 1)
        ), proj4string=crsLL
    )
    adjust = mapply(
        function(crs){
          pn2 = spTransform(
              pointNorth,
              crs
              )
          pnDist =apply(pn2@coords,2,diff)
          -atan(pnDist[1]/pnDist[2])*360/(2*pi)
        },
        crs=rotatedCRS
    )

    adjust[adjust<0] =
        360+adjust[adjust<0]
    
    rotatedProj4stringAdj = 
        omercProj4string(
            lon=theCentre[1],
            lat=theCentre[2], 
            angle=angle,
            inverseAngle=adjust
    )
    rotatedCrsAdj = lapply(rotatedProj4stringAdj, CRS)
  } else {
    rotatedCrsAdj = rotatedCRS
  }
  
  
  if(is.numeric(x) | !haveRgdal) {
      if(length(rotatedCrsAdj)==1)
        rotatedCrsAdj = rotatedCrsAdj[[1]]
      return(rotatedCrsAdj)
  }

  # find the optimal rotatino
  # for a small bounding box
    xTrans = mapply(
        function(CRSobj) {
          abs(prod(apply(bbox(
               spTransform(x, CRSobj)           
                          ), 1, diff)
              ))
        },
        CRSobj=rotatedCRS
        )
    rotatedCrsAdj = rotatedCrsAdj[[
        which.min(xTrans)
        ]]
        
        
  result = spTransform(x, rotatedCrsAdj)
        
  result
  
}
  
  