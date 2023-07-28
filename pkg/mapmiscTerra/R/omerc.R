omercProj4string = function(
    lon, lat, angle, 
    x=0,y=0, inverseAngle=0,
    scale=1,
    ellps='WGS84', units='m',
    datum='WGS84',
    crs=TRUE) {
  
#  negAngle = angle<0
#  angle[negAngle] = 360 + angle[negAngle]
#  angle[angle==90]=89

  x = rep_len(x, length(angle))
  y = rep_len(y, length(angle))
  scale = rep_len(scale, length(angle))
  lat = rep_len(lat, length(angle))
  lon = rep_len(lon, length(angle))

  
  whichZeros = (angle==0) 
  which90 = abs(angle)==90
  
  
  result = paste(
      "+proj=omerc",
      " +lat_0=", lat,
      " +lonc=", lon,
      " +alpha=", angle, 
      " +k_0=", scale, 
      " +x_0=", x,
      " +y_0=", y,
      " +gamma=", inverseAngle,
      " +ellps=", ellps,
      " +units=", units,
      ' +datum=', datum,
      sep="")
  
  if(any(whichZeros)) {
    result[whichZeros] = paste(
        "+proj=tmerc",
        " +lat_0=", lat[whichZeros] ,
        " +lon_0=", lon[whichZeros] ,
        " +k=", scale[whichZeros] , 
        " +x_0=", x[whichZeros] ,
        " +y_0=", y[whichZeros] ,
        " +ellps=", ellps,
        " +units=", units,
        ' +datum=', datum,
        sep="")
  }

  
  if(crs) {
    result = lapply(result, crs)
  }
  
  result
}

omerc = function(
    x, angle=0, 
    post=c('none', 'north', 'wide','tall'),
    preserve=NULL,
    ellipse=TRUE
) {
  

  digits=3
  angleOrig = angle
  post = post[1]
  
  if(is.numeric(post)){
    inverseAngle = rep_len(post, length(angle))
    post='none'
  } else {
    inverseAngle=rep(0, length(angle))
    post=post[seq(1,len=length(post))]
  }
  scale = rep(1, length(angle))
  objectiveResult = NULL
  
  # convert all angles to -90 to 90 range
  southAngle = ( abs(angle)>90) & (abs(angle) < 270 )
  
  
  angle[southAngle] = angle[southAngle] + 180
  inverseAngle[southAngle] = inverseAngle[southAngle] + 180
  
  angleC = exp(1i*2*pi*(angle/360))
  angle = round(Arg(angleC)*360/(2*pi), digits)
  
#  the90s = abs(angle)==90
#  inverseAngle[the90s] = inverseAngle[the90s] + angle[the90s]
#  angle[the90s] = 0

  # create the centre point
  if(!is.numeric(x)){
    # centre of the bounding box
    theCentre = terra::ext(x)
    theCentre = c(terra::xmin(theCentre), terra::ymin(theCentre)) +
        diff(as.vector(theCentre))[-2]/2
  } else {
    theCentre = x[1:2]
    x = vect(matrix(x, ncol=2), crs=crsLL)
    # use preserve for finding best bounding box
    if(!is.null(preserve)){
      x = preserve
    } else {
      # otherwise can't find best bounding box
      if(post %in% c('wide', 'tall'))
        post = 'none'
    }
  }
  theCentre = round(theCentre, digits)
  
  crs = crs(x)
  if(is.na(crs)){
    crs = crsLL
  }

  # convert the centre to LL if necessary
  if(!terra::is.lonlat(crs)) {
      theCentre = vect(
          t(theCentre[1:2]),
          crs=crs
          )
      theCentre = project(theCentre, crsLL)
      theCentre = as.vector(terra::crds(theCentre))
  } # crs not LL

  
  # projections without inverse rotation
  rotatedCRS = 
        omercProj4string(
            lon=theCentre[1],
            lat=theCentre[2], 
            angle=angle)
    
    
    # make sure theCentre is at the origin
   theCentreSp = vect(t(theCentre), crs=crsLL)
   newxy = simplify2array(lapply(rotatedCRS, function(qq){
          drop(terra::crds(project(theCentreSp, qq)))
        }))
    newxy = round(newxy, digits)
    
    rotatedCRS = omercProj4string(
        lon=theCentre[1],
        lat=theCentre[2], 
        angle=angle,
        x=-newxy[1,],
        y=-newxy[2,])


    
    # preserve distances between points
    if(!is.null(preserve)) {
      # convert to LL
      if(!terra::is.lonlat(crs(preserve))){
        preserve = project(preserve, crsLL)
      }
      # great circle distance
	      distGS = as.matrix(terra::distance(preserve)*1000)
      theLower = lower.tri(distGS, diag=FALSE)
      distGS = distGS[theLower]
       # euclidean distance for each projection
      distEu = unlist(lapply(rotatedCRS,
              function(crs){
                mean(
                    as.matrix(terra::distance(
                        project(preserve, crs)
                ))[theLower]/distGS,
                    na.rm=TRUE)
              }
          ))
      # add scaling to CRS
      rotatedCRS = 
          omercProj4string(
              lon=theCentre[1],
              lat=theCentre[2], 
              angle=angle,
              x=-newxy[1,],
              y=-newxy[2,],
              scale=round(1/distEu, digits))
      
      # find ratio of Euclidean to great circle distances
      distEuSsq = unlist(
          lapply(rotatedCRS,
              function(crs){
                sqrt(mean(
                        (terra::distance(project(preserve, crs)
                              )/distGS
                              -1)^2,
                        na.rm=TRUE))
              }
          )
      )

      minDist= which.min(distEuSsq)
      objectiveResult=list(
          x = angle,
          y = distEuSsq
      )
      
      angle=angle[minDist]
      scale=round(1/distEu[minDist], digits)
      inverseAngle = inverseAngle[minDist]
      newxy = newxy[,minDist,drop=FALSE]
      
      
      
     } # end preserve
    
    
    # find the optimal rotation
    # for a small bounding box
    # if x is not numeric, and more than one crs, not preserving distances
    if(!is.numeric(x) & (length(rotatedCRS)>1) & is.null(preserve)){
      
      # assume the best 'tall' bounding box will be found
      if(post=='tall')
        post = 'none'
      if(post=='wide'){
        post='none'
        inverseAngle = rep(90, length(angle))
      }
      
      xTrans = mapply(
          function(CRSobj) {
            prod(abs(diff(terra::ext(
                            project(x, CRSobj)           
                        ))[-2]
                ))
          },
          CRSobj=rotatedCRS
      )
      

      objectiveResult=list(
          x = angle,
          y = xTrans)

      minX = which.min(xTrans)
      angle = angle[minX]
      scale = scale[minX]
      inverseAngle = inverseAngle[minX]
      newxy = newxy[,minX,drop=FALSE]
      


      
    } # end smallest bounding box

    rotatedCRS = omercProj4string(
        lon=theCentre[1],
        lat=theCentre[2], 
        angle=angle,
        scale=scale,
        x=-newxy[1,],
        y=-newxy[2,],
        inverseAngle=inverseAngle
    )
    
    if(post=='none') {
      if(length(rotatedCRS)==1) {
        rotatedCRS = rotatedCRS[[1]]
      }
      if(!is.null(objectiveResult))
        attributes(rotatedCRS)$obj = objectiveResult
#      return(rotatedCRS)
    }
    
# otherwise find a better inverse angle



    
  # find an inverse rotation to preserve north
  if(post=='north') {
    # ignore any existing inverseAngle
    rotatedCRS = omercProj4string(
        lon=theCentre[1],
        lat=theCentre[2], 
        angle=angle,
        scale=scale,
        x=-newxy[1,],
        y=-newxy[2,]
    )
    # a pair of points which should be 
  # north-south of each other
    pointNorth = vect(
        rbind(
            theCentre,
            theCentre + c(0, 0.1)
        ), crs=crsLL
    )
 
    adjust = mapply(
        function(crs){
          pn2 = project(
              pointNorth,
              crs
              )
          pn2 = vect(round(terra::crds(pn2), digits), crs=crs)
              
          # find their distance in new projection
          pnDist =apply(terra::crds(pn2),2,diff)
          # and the angle they are from North-South
          -atan(pnDist[1]/pnDist[2])*360/(2*pi)
        },
        crs=rotatedCRS
    )
 
    inverseAngle = round(adjust, digits)
 
  } # end post n orth
  
  if(post=='wide' | post=='tall' & (length(angle)==1)) {
    Sgamma = seq(-90,90,)
    # ignore any existing inverseAngle
    rotatedCRS = omercProj4string(
        lon=theCentre[1],
        lat=theCentre[2], 
        angle=angle,
        scale=scale,
        x=-newxy[1,],
        y=-newxy[2,],
        inverseAngle = Sgamma
    )
    bbarea = mapply(
        function(CRSobj) {
          thebb = diff(as.vector(terra::ext(
              project(x, CRSobj)           
          )))[-2]
          c(area= abs(prod(thebb)),
              ratio = as.numeric(abs(thebb[2]/thebb[1]))
          )
        },
        CRSobj=rotatedCRS
    )
    bbarea = rbind(bbarea, inverseAngle=Sgamma)
    if(post=='wide'){
      bbarea = bbarea[,bbarea['ratio',]<=1]
    } else {
      bbarea = bbarea[,bbarea['ratio',]>=1]
    }
    inverseAngle = bbarea[
        'inverseAngle',
        which.min(bbarea['area',])]
  }
  
  
  # create new proj4string
# with xy offset and inverse rotation
  rotatedCRS = 
      omercProj4string(
          lon=theCentre[1],
          lat=theCentre[2], 
          angle=angle,
          inverseAngle=inverseAngle,
          x=-newxy[1,],
          y=-newxy[2,],
          scale=scale)
  
  
  if(length(rotatedCRS)==1) {
    rotatedCRS = rotatedCRS[[1]]
  }

  if(!is.null(objectiveResult)) {
    attributes(rotatedCRS)$obj = objectiveResult
  }
  

  equatorMeetsGreatCircle = geosphere::gcIntersectBearing(theCentre, angle, c(0,0), 90)
  if(any(is.na(equatorMeetsGreatCircle)))
  equatorMeetsGreatCircle = geosphere::gcIntersectBearing(theCentre, angle, c(0,0), -90)
  if(any(is.na(equatorMeetsGreatCircle))) {
    theY = 0
  } else {
  equatorMeetsGreatCircle = vect(matrix(equatorMeetsGreatCircle, ncol=2, byrow=TRUE), crs=crsLL)
  equatorMeetsGreatCircleProj = project(equatorMeetsGreatCircle, rotatedCRS)
  theY = terra::crds(equatorMeetsGreatCircleProj)[,2]
  theY = theY[which.min(abs(theY))]
  }


#  theEllipse = terra::as.polygons(ext(c(-1.6e7, 1.6e7, -2e7, 2e7) + theY*c(0,0,1,1)), rotatedCRS)
  theEllipse = terra::as.polygons(terra::ext(extentMerc + theY*c(0,0,1,1)), rotatedCRS)


  if(ellipse) {
  theBox = llCropBox(crs=rotatedCRS, ellipse = theEllipse, 
      buffer.width=40*1000, densify.interval = 20*1000, 
        crop.distance = Inf, crop.poles = FALSE, crop.leftright=TRUE,
      remove.holes=TRUE)
} else {
  theBox = list(ellipse=theEllipse)
}
#  plot(project(worldMap, crsLL));plot(theBox$crop, col='red', add=TRUE)   

  attributes(rotatedCRS)[names(theBox)] = theBox


  return(rotatedCRS)
}
