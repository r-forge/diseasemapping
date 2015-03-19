omercProj4string = function(
    lon, lat, angle, 
    x=0,y=0, inverseAngle=0,
    scale=1,
    ellps='WGS84', units='m',
    crs=TRUE) {
  
#  negAngle = angle<0
#  angle[negAngle] = 360 + angle[negAngle]
#  angle[angle==90]=89

  whichZeros = angle==0
  
  result = paste(
      "+proj=omerc",
      " +lat_0=", lat,
      " +lonc=", lon,
      " +alpha=", angle, 
      " +k=", scale, 
      " +x_0=", x,
      " +y_0=", y,
      " +gamma=", inverseAngle,
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
  if(crs) 
    result = lapply(result, CRS)
  
  result
}

omerc = function(
    x, angle=0, undo=FALSE,
    preserve=NULL
) {
  
  digits=3 # for rounding coordinates
  angle = round(angle, digits)
  
  angle = angle[! (angle %in% (90*c(-1,1,2)))]
  if(!length(angle)){
    warning('angle cant be -90, 90 or 120')
  }
  
  
  crs = projection(x)
  if(is.na(crs)){
    crs = crsLL
  }
  
  if(is.character(crs))
    crs = CRS(crs)
  
  haveRgdal = requireNamespace('rgdal', quietly=TRUE)

  # create the centre point
  if(!is.numeric(x)){
    # centre of the bounding box
    theCentre = bbox(x)
    theCentre = theCentre[,'min'] +
      apply(theCentre, 1, diff)/2
  } else {
    theCentre = x[1:2]
  }
  theCentre = round(theCentre, digits)

  # convert the centre to LL if necessary
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

  # create strings for projections
  rotatedCRS = 
      omercProj4string(
      lon=theCentre[1],
      lat=theCentre[2], 
      angle=angle)
  
  # some refinements if gdal is available
  if(haveRgdal) {
    # make sure theCentre is at the origin
   theCentreSp = SpatialPoints(t(theCentre), proj4string=crsLL)
   newxy = simplify2array(lapply(rotatedCRS, function(qq){
          drop(spTransform(theCentreSp, qq)@coords)
        }))
    newxy = round(newxy)
    rotatedCRS = omercProj4string(
        lon=theCentre[1],
        lat=theCentre[2], 
        angle=angle,
        x=-newxy[1,],
        y=-newxy[2,])
  } else{
    # no gdal, don't re-centre
    newxy = matrix(0,2,1)
  }
  
  # find an inverse rotation to preserve north
  if(undo & haveRgdal) {
    # a pair of points which should be 
  # north-south of each other
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
              
          # find their distance in new projection
          pnDist =apply(pn2@coords,2,diff)
          # and the angle they are from North-South
          -atan(pnDist[1]/pnDist[2])*360/(2*pi)
        },
        crs=rotatedCRS
    )
    adjust = round(adjust, digits)
    #    adjust[adjust<0] =
#        360+adjust[adjust<0]
    
    # make sure the centre is still correct
    # redo the adjustment

# create new proj4strings 
# with inverse rotation
# but no centering
  rotatedCrsAdj = 
        omercProj4string(
            lon=theCentre[1],
            lat=theCentre[2], 
            angle=angle,
            inverseAngle=adjust,
            x=0,
            y=0
        )
    
    # find coordinates of centre point
      newxy = simplify2array(lapply(rotatedCrsAdj, function(qq){
                drop(spTransform(theCentreSp, qq)@coords)
              }))
      newxy = round(newxy)
      
    # create new proj4string
# with xy offset and inverse rotation
  rotatedCrsAdj = 
          omercProj4string(
              lon=theCentre[1],
              lat=theCentre[2], 
              angle=angle,
              inverseAngle=adjust,
              x=-newxy[1,],
              y=-newxy[2,])
   } else {
     # no gdal or no undo, 
  # set adjusted CRS to unadjusted CRS
    rotatedCrsAdj = rotatedCRS
  }
  
  # preserve distances between points
  if(haveRgdal & !is.null(preserve)) {
    # convert to LL
    if(!isLonLat(projection(preserve))){
      preserve = spTransform(preserve, crsLL)
    }
    # great circle distance
    distGS = spDists(preserve, longlat=TRUE)*1000
    theLower = lower.tri(distGS, diag=FALSE)
    distGS = distGS[theLower]
    # euclidean distance for each projection
    distEu = unlist(lapply(rotatedCrsAdj,
        function(crs){
          mean(
              spDists(spTransform(preserve, crs), 
                  longlat=FALSE)[theLower]/distGS,
              na.rm=TRUE)
        }
    ))
    # add scaling to CRS
  rotatedCrsAdj = 
      omercProj4string(
          lon=theCentre[1],
          lat=theCentre[2], 
          angle=angle,
          inverseAngle=adjust,
          x=-newxy[1,],
          y=-newxy[2,],
          scale=round(1/distEu, digits))

    # find ratio of Euclidean to great circle distances
    distEu = unlist(
      lapply(rotatedCrsAdj,
      function(crs){
        sqrt(mean(
            (spDists(spTransform(preserve, crs), 
                longlat=FALSE)[theLower]/distGS
            -1)^2,
            na.rm=TRUE))
        }
      )
    )
    # plot(angle, distEu, log='y')
  # select CRS with best preserved distances
  rotatedCrsAdj = rotatedCrsAdj[
      which.min(distEu)
  ]
  attributes(rotatedCrsAdj[[1]])$obj=list(
      x = angle,
      y = distEu
      )
  
  } # end preserve
  
  # find the optimal rotatino
  # for a small bounding box
  # if x is not numeric, have gdal, and more than one crs
  if(!is.numeric(x) & haveRgdal & (length(rotatedCrsAdj)>1)){
    xTrans = mapply(
      function(CRSobj) {
        abs(prod(apply(bbox(
                        spTransform(x, CRSobj)           
                    ), 1, diff)
            ))
      },
      CRSobj=rotatedCRS
    )
    rotatedCrsAdj = rotatedCrsAdj[
      which.min(xTrans)
    ]
    attributes(rotatedCrsAdj[[1]])$obj=list(
        x = angle,
        y = xTrans
    )

  } # end smallest bounding box

  if(length(rotatedCrsAdj)==1)
     rotatedCrsAdj = rotatedCrsAdj[[1]]

  rotatedCrsAdj
}
