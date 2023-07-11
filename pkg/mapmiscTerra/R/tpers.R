# http://proj4.org/projections/tpers.html
tpers = function(x,
  hKm = 100*1000, tilt = -10,
  offset=c(0,0), 
  axis='enu') {
  
  if(!terra::is.lonlat(x)) x = project(x, crsLL)

  myCrs = crs(paste(
      "+proj=tpers +h=",
      hKm*1000,
      " +lat_0=",
      x@coords[1,2],
      " +lon_0=",
      x@coords[1,1],
      " +azi=",
      geosphere::bearing(x@coords[1,],x@coords[2,])[1],
      " +tilt=", tilt,
      " +ellps=WGS84 +axis=", axis,
      " +x_0=", offset[1],
      " +y_0=", offset[2],
      sep=""))
  
  cropBox = llCropBox(crs=myCrs)
  
  attributes(myCrs)$ellipse = cropBox$ellipse
  attributes(myCrs)$crop = cropBox$crop
   attributes(myCrs)$regionLL = cropBox$regionLL

  
  myCrs
}
