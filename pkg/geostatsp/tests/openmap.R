# only do the following if running unix (because nsl is available)
# and if the OpenStreetMap.org web site can be accessed

myraster = raster(matrix(0,10,10),xmn=8,xmx=18,ymn=0,ymx=10, crs="+proj=longlat")
values(myraster) = seq(0,1,len=ncell(myraster))
myPoints = SpatialPoints(myraster, proj4string=CRS(proj4string(myraster)))[
		seq(1,ncell(myraster),len=5)]

plot(myraster)
points(myPoints)

utmproj = "+proj=utm +zone=32" 
myrasterUTM = projectRaster(myraster, crs=utmproj)
myPointsUTM = spTransform(myPoints, CRS(utmproj))
plot(myrasterUTM)
points(myPointsUTM)

myPointsMercator = spTransform(myPoints, CRS("+proj=merc +a=6378137 +b=6378137 +lat_ts=0.0 +lon_0=0.0 +x_0=0.0 +y_0=0
						+k=1.0 +units=m"))


if(.Platform$OS.type=="unix") {
	if(length(utils::nsl("www.OpenStreetMap.org"))) {
		par(mar=c(0,0,0,0))
		mytiles = openmap(myraster, minNumTiles=2,  type="osm")
		plot(myraster)
		plot(mytiles, add=TRUE)
		points(myPoints,col='red')
		
		mytiles = openmap(extent(myraster), minNumTiles=2,  type="osm")
		plot(myPointsMercator)
		plot(mytiles, add=TRUE)
		points(myPointsMercator,col='red')
		
		mytiles = openmap(bbox(myraster), minNumTiles=2,   type="osm")
		plot(myPointsMercator)
		plot(mytiles, add=TRUE)
		points(myPointsMercator,col='red')
		
		mytiles = openmap(myPoints, minNumTiles=2,  type="osm")
		plot(myPoints)
		plot(mytiles, add=TRUE)
		points(myPoints,col='red')
		
		
		mytiles = openmap(myrasterUTM, minNumTiles=2,  type="osm")
		plot(myrasterUTM)
		plot(mytiles, add=TRUE)
		points(myPointsUTM,col='red')
		
		
		mytiles = openmap(myPointsUTM, minNumTiles=2,  type="osm")
		plot(myPointsUTM)
		plot(mytiles, add=TRUE)
		points(myPointsUTM,col='red')
		
	}
}