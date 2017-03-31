
library('mapmisc')
library('rgdal')

okinawa = mapmisc::geocode("Okinawa")
y = mapmisc::geocode("Hong Kong")
x = mapmisc::geocode("Guam")
korea = raster::getData(
  "GADM", country='KOR', level=0
)


myProj = #mapmisc:::
tpers(x,y,hKm = 6*1000,tilt =-25)

data("wrld_simpl", package='maptools')

wrld4 = wrapPoly(wrld_simpl, myProj)

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
#  path='nrcan')
#path='https://sat01.maps.yandex.net/tiles?l=sat&v=1.35.0&')
path='https://services.arcgisonline.com/ArcGIS/rest/services/Ocean/World_Ocean_Base/MapServer/tile/')

par(bg='black')
map.new(myMap, buffer=-c(0,0,3,0)*1000*1000)
plotRGB(myMap, add=TRUE, colNA='black')

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
plot(koreaT, add=TRUE, border='yellow')  



