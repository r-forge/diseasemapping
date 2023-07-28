# http://proj4.org/projections/tpers.html
tpers = function(x,
  hKm = 100*1000, tilt = -10,
  azi,
  offset=c(0,0), 
  axis='enu') {

  if(is.vector(x)) x = vect(matrix(x, ncol=2), crs=crsLL)
  
  if(!terra::is.lonlat(x)) x = project(x, crsLL)


  if(missing(azi)) {
  if(length(x)==1) {
    azi = 0
  } else {
    azi = geosphere::bearing(terra::crds(x)[1,],terra::crds(x)[2,])[1]    
  }
  }

  myCrs = terra::crs(paste(
      "+proj=tpers +h=",
      hKm*1000,
      " +lat_0=",
      terra::crds(x)[1,2],
      " +lon_0=",
      terra::crds(x)[1,1],
      " +azi=", azi,
      " +tilt=", tilt,
      " +ellps=WGS84 +axis=", axis,
      " +x_0=", offset[1],
      " +y_0=", offset[2],
      sep=""))
  
  cropBox = llCropBox(crs=myCrs, crop.poles=FALSE, remove.holes=FALSE,  buffer.width=10*1000, densify.interval = 2000)
  
  attributes(myCrs)$ellipse = cropBox$ellipse
  attributes(myCrs)$crop = cropBox$crop
   attributes(myCrs)$regionLL = terra::fillHoles(cropBox$crop)

  
  myCrs
}
