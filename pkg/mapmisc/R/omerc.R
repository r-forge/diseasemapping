omercProj4string = function(
    lon, lat, angle, 
    x=0,y=0, inverseAngle=0,
    scale=1,
    ellps='WGS84', units='m') {
  
#  negAngle = angle<0
#  angle[negAngle] = 360 + angle[negAngle]
#  angle[angle==90]=89

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
    x, angle=0, undo=FALSE,
    preserve=NULL
) {
  
  digits=3
  angle = round(angle, digits)
  
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
  theCentre = round(theCentre, digits)
   
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
  
  # make sure theCentre is at the origin
  if(haveRgdal) {
   theCentreSp = SpatialPoints(t(theCentre), proj4string=crsLL)
   newxy = simplify2array(lapply(rotatedCRS, function(qq){
          drop(spTransform(theCentreSp, qq)@coords)
        }))
    newxy = round(newxy)
    rotatedProj4string = 
    omercProj4string(
        lon=theCentre[1],
        lat=theCentre[2], 
        angle=angle,
        x=-newxy[1,],
        y=-newxy[2,])
    rotatedCRS = lapply(rotatedProj4string, CRS)
  } else{
    newxy = matrix(0,2,1)
  }
  if(undo & haveRgdal) {
    pointNorth = SpatialPoints(
        rbind(
            theCentre,
            theCentre + c(0, 0.1)
        ), proj4string=crsLL
    )
    adjust = mapply(
        function(crs){
          pn2 = spTransform(
              pointNorth,
              crs
              )
          pn2@coords = round(pn2@coords)
              
          pnDist =apply(pn2@coords,2,diff)
          -atan(pnDist[1]/pnDist[2])*360/(2*pi)
        },
        crs=rotatedCRS
    )

#    adjust[adjust<0] =
#        360+adjust[adjust<0]
    
    rotatedProj4stringAdj = 
        omercProj4string(
            lon=theCentre[1],
            lat=theCentre[2], 
            angle=angle,
            inverseAngle=adjust,
            x=-newxy[1,],
            y=-newxy[2,]
    )
    rotatedCrsAdj = lapply(rotatedProj4stringAdj, CRS)
    
      # redo the adjustment
    rotatedProj4stringAdj = 
        omercProj4string(
            lon=theCentre[1],
            lat=theCentre[2], 
            angle=angle,
            inverseAngle=adjust,
            x=0,
            y=0
        )
    rotatedCrsAdj = lapply(rotatedProj4stringAdj, CRS)
    
      newxy = simplify2array(lapply(rotatedCrsAdj, function(qq){
                drop(spTransform(theCentreSp, qq)@coords)
              }))
      newxy = round(newxy)
      rotatedProj4stringAdj = 
          omercProj4string(
              lon=theCentre[1],
              lat=theCentre[2], 
              angle=angle,
              inverseAngle=adjust,
              x=-newxy[1,],
              y=-newxy[2,])
      rotatedCrsAdj = lapply(rotatedProj4stringAdj, CRS)
   } else {
    rotatedCrsAdj = rotatedCRS
  }
  
  if(haveRgdal & !is.null(preserve)) {
    if(!isLonLat(projection(preserve))){
      preserve = spTransform(preserve, crsLL)
    }
    distGS = spDists(preserve, longlat=TRUE)*1000
    
    distEu = lapply(rotatedCRS,
        function(crs){
          mean(
              spDists(spTransform(preserve, crs), longlat=FALSE)/distGS,
              na.rm=TRUE)
        }
    )
  rotatedCrsAdj = mapply( 
      function(crs,scale){
        CRS(gsub(
                "\\+k=1", 
                paste("+k=", round(1/scale,digits),sep=''),
                as.character(crs)))
        
      },
      crs=rotatedProj4stringAdj,
      scale=distEu)
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
        

        rotatedCrsAdj

}
  
  