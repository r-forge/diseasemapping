
library('mapmisc')
library('rgdal')

options(mapmiscCachePath = system.file('extdata', package='mapmisc'))
options(mapmiscVerbose=TRUE)

okinawa = mapmisc::geocode("Okinawa")
hk = y = mapmisc::geocode("Hong Kong")
guam = mapmisc::geocode("Guam")
x = SpatialPoints(cbind(130,15), proj4string=crsLL)


myProj = tpers(x=bind(x,y), hKm = 2*1000,
  tilt =-15, axis='seu')

if(length(attributes(myProj)$ellipseSafeLL)>2) {
  
myMap1 = openmap(
  x=attributes(myProj)$regionLL[1,],
  crs=crsMerc,
  verbose=TRUE,
  zoom=3,
path='https://sat01.maps.yandex.net/tiles?l=sat&v=1.35.0&')
#  path='https://services.arcgisonline.com/ArcGIS/rest/services/Ocean/World_Ocean_Base/MapServer/tile/')

myMap2 = openmap(
  x=attributes(myProj)$regionLL[2,],
  zoom=attributes(myMap1)$tiles$zoom,
  crs=crsMerc,
  verbose=TRUE, zoom=3,
path='https://sat01.maps.yandex.net/tiles?l=sat&v=1.35.0&')
#  path='https://services.arcgisonline.com/ArcGIS/rest/services/Ocean/World_Ocean_Base/MapServer/tile/')

mapRaster = raster(
  extent(attributes(myProj)$ellipse), 
  res = mean(res(myMap1)/1.5),
  crs=myProj
)

myMap1a = projectRaster(
  from=myMap1, to=mapRaster, method='ngb')
myMap2a = projectRaster(
  from=myMap2, to=mapRaster, method='ngb')

myMap = mosaic(myMap1a, myMap2a, fun=mean)
} else {
  myMap = openmap(
    x=attributes(myProj)$regionLL,
    buffer=1, zoom=3,
    crs=myProj,
    verbose=TRUE,
    path='https://sat01.maps.yandex.net/tiles?l=sat&v=1.35.0&')
#    path='https://services.arcgisonline.com/ArcGIS/rest/services/Ocean/World_Ocean_Base/MapServer/tile/')
  }

data("wrld_simpl", package='maptools')

wrld4 = wrapPoly(x=wrld_simpl, crs=myProj)



xProj = spTransform(
  guam, myProj
)
yProj = spTransform(hk, myProj)
okinawaProj = spTransform(
  okinawa, myProj
)





map.new(myProj, bg='black', buffer=2*100*1000)
plotRGB(myMap, add=TRUE, bgalpha=0)

points(xProj, col='red', pch=16, cex=2)
text(xProj@coords, 
  as.character(xProj$originalPlace),
  pos=2, col='yellow')

points(okinawaProj, col='red', pch=16, cex=2)
text(okinawaProj@coords, 
  as.character(okinawaProj$originalPlace),
  pos=2, col='yellow')

points(yProj, col='red', pch=16, cex=2)
text(yProj@coords, 
  as.character(yProj$originalPlace),
  pos=2, col='yellow')

plot(wrld4, add=TRUE)

gridlinesWrap(myProj, col='yellow', 
  north=seq(-40,40,by=10), easts=seq(-150,180,by=10),
  lty=2)



