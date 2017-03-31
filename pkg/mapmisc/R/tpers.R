# http://proj4.org/projections/tpers.html
tpers = function(x,y,
  hKm = 100*1000, tilt = -10) {
  myCrs = CRS(paste(
      "+proj=tpers +h=",
      hKm*1000,
      " +lat_0=",
      x@coords[1,'latitude'],
      " +lon_0=",
      x@coords[1,'longitude'],
      " +azi=",
      geosphere::bearing(x,y)[1],
      " +tilt=", tilt,
      " +ellps=WGS84", 
      sep=""))
  
  attributes(myCrs)$retain =
    llCropBox(myCrs, res=1)$crop
  
  polarLimit = 91
  extentLL = as(
    extent(-180,180,-polarLimit,polarLimit),
    'SpatialPolygons')
  extentLL@proj4string = crsLL
  
  attributes(myCrs)$crop =
    suppressWarnings(rgeos::gBuffer(
        rgeos::gDifference(extentLL, 
          attributes(myCrs)$retain, byid=FALSE),
        width=0.5
      ))
  attributes(myCrs)$retain =
    suppressWarnings(
      rgeos::gDifference(attributes(myCrs)$retain, 
        attributes(myCrs)$crop))
  
  attributes(myCrs)$ellipse = 
    rgeos::gBuffer(spTransform(
    attributes(myCrs)$retain, myCrs
  ), byid=TRUE, width=70*1000)

  attributes(myCrs)$ellipse = 
      rgeos::gUnaryUnion(attributes(myCrs)$ellipse)

  myCrs
}
