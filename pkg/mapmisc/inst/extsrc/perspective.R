
library('mapmisc')
library('rgdal')

okinawa = mapmisc::geocode("Okinawa")
y = mapmisc::geocode("Hong Kong")
x = mapmisc::geocode("Guam")
korea = raster::getData(
  "GADM", country='KOR', level=0
)


# http://proj4.org/projections/tpers.html
(myProj = CRS(paste(
        "+proj=tpers +h=",
        6*1000*1000,
        " +lat_0=",
        x@coords[1,'latitude'],
        " +lon_0=",
        x@coords[1,'longitude'],
        " +azi=",
        geosphere::bearing(x,y),
        " +tilt=-25",
        " +ellps=WGS84", 
        sep="")))

koreaT = spTransform(korea,
  myProj)  
xProj = spTransform(
  x, myProj
)
yProj = spTransform(y, myProj)
okinawaProj = spTransform(
  okinawa, myProj
)



myMap = openmap(raster( 
  extent(100,180,-10,70), crs=crsLL),
  fact=1.5, maxTiles=24,
  crs = myProj,
  verbose=TRUE,
#  path=mPath)
#  path = 'https://a.tile.hosted.thunderforest.com/komoot-2/')
path='https://sat01.maps.yandex.net/tiles?l=sat&v=1.35.0&')
#  path='https://server.arcgisonline.com/ArcGIS/rest/services/World_Physical_Map/MapServer/tile/')
#path='https://server.arcgisonline.com/ArcGIS/rest/services/World_Imagery/MapServer/tile/')
#path='https://services.arcgisonline.com/ArcGIS/rest/services/Ocean/World_Ocean_Base/MapServer/tile/')

map.new(myMap, 
  buffer=-c(0,0,3,0)*1000*1000)
plotRGB(myMap, add=TRUE)

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

plot(koreaT, add=TRUE, border='yellow')  



