# crsLL = CRS("+epsg:4326")
crsLL = CRS("+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0") 

#  crsMerc =CRS("+init=epsg:3857") # mercator projection
crsMerc = CRS("+proj=merc +a=6378137 +b=6378137 +lat_ts=0.0 +lon_0=0.0 +x_0=0.0 +y_0=0 +k=1.0 +ellps=WGS84 +units=m +nadgrids=@null +no_defs")

crsMercSphere = CRS("+proj=merc +ellps=sphere +units=m")

crsSphere = CRS("+proj=longlat +ellps=sphere")

llLim = atan(sinh(pi))*360/(2*pi)
openmapExtentLL = extent(-180, 180,-llLim,llLim)

if(FALSE) {
openmapExtentMercSphere = extent(projectExtent(raster(openmapExtentLL, crs=crsSphere), crsMercSphere))
openmapExtentMercSphere = extent(round(as.vector(openmapExtentMercSphere)))
openmapExtentMercSphere = as.vector(openmapExtentMercSphere)
dump('openmapExtentMercSphere', file='')

openmapExtentMerc = extent(projectExtent(raster(openmapExtentLL, crs=crsSphere), crsMerc))
openmapExtentMerc = extent(round(as.vector(openmapExtentMerc)))
openmapExtentMerc = as.vector(openmapExtentMerc)
dump('openmapExtentMerc', file='')

}
openmapExtentMercSphere <-
    extent(c(-20015077, 20015077, -20015077, 20015077))

openmapExtentMerc <-
   extent( c(-20037508, 20037508, -19994875, 19994875))

.getExtent = function(x, crs=NA, extend=0, crsOut = crsMercSphere) {
  
  # find the CRS's
  crsUse = projection(x)
  
#	if(is.logical(crs)) { # it's probably NA,
  if(all(is.na(crs))) {
    crs=crsLL
  }
#	if(is.logical(crsUse)) { # it's probably NA,
  if(all(is.na(crsUse))) {
    crsUse=crs
  }
  
  if(identical(crsUse ,"NA"))
    crsUse = crs
  
  # if x is numeric, transform to extent  
  tileEps = sqrt(.Machine$double.neg.eps)
  if(is.numeric(x))
    if(length(x)==2)
      x = extent(x[1]-tileEps, xmax=x[1]+tileEps, ymin=x[2]-tileEps, ymax=x[2]+tileEps)
  if(requireNamespace('rgdal', quietly=TRUE)) {
    x = as(
        extent(x), 
        'SpatialPoints'
    )
    crs(x) = crsUse
    result =  raster::extend(extent(spTransform(x, crsOut)), extend)
  } else {
    # no rgdal, try to use raster, doesn't always work
    x = raster(extent(x), nrows=100,ncols=100, crs=crsUse)
    x = raster::extend(x, extend(extent(x), extend))
    result = extent(projectExtent(x, crsOut))
  }
  
  result
}

cropExtent = function(x,y){
  xmin(x) = pmax(xmin(x), xmin(y))
  xmax(x) = pmin(xmax(x), xmax(y))
  ymin(x) = pmax(ymin(x), ymin(y))
  ymax(x) = pmin(ymax(x), ymax(y))
  x
}

