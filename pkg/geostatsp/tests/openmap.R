
library(geostatsp)
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


myplot = function(first,second=first) {
	pdf(tempfile("osmplot", tmpdir=".", fileext=".pdf"))
	par(mar=c(0,0,0,0))
	plot(first)
	plot(mytiles, add=TRUE)
	plot(second,add=TRUE,col='blue')
	points(mycities,col='red')
	text(mycities, labels=mycities$name, col='red',pos=4)
	dev.off()
}

# only do the following if running unix (because nsl is available)
# and if the OpenStreetMap.org web site can be accessed
if(exists("nsl", where="package:utils")) {
	if(length(utils::nsl("www.OpenStreetMap.org"))& 
			length(utils::nsl("ws.geonames.org"))) {


		# raster, result will be in project of the raster (long-lat)
		mytiles = openmap(extend(myraster,1), minNumTiles=2)
		mycities = GNcities(extend(myraster,1),max=5)
		myplot(myraster, myPoints)
		
		# extent, tiles will be mercator
		mytiles = openmap(extent(myraster), minNumTiles=2)
		# cities will be long=lat
		mycities = GNcities(extent(myraster),max=5,lang="fr")
		# so change to mercator
		mycities = spTransform(mycities, CRS(proj4string(myPointsMercator)))
		myplot(myPointsMercator)
		
		# give the bbox, again tiles mercator, cities long lat
		mytiles = openmap(bbox(myraster), minNumTiles=2)
		mycities = GNcities(bbox(myraster),max=5)
		mycities = spTransform(mycities, CRS(proj4string(myPointsMercator)))
		myplot(myPointsMercator)
		
		# give points, result is CRS of points (long-lat)
		mytiles = openmap(myPoints, minNumTiles=2)
		mycities = GNcities(myPoints,max=5,lang="es")
		myplot(myPoints)
		
		# UTM raster
		mytiles = openmap(myrasterUTM, minNumTiles=2)
		mycities = GNcities(myrasterUTM,max=5)
		myplot(myrasterUTM, myPointsUTM)
		
		# utm points
		mytiles = openmap(myPointsUTM, minNumTiles=2)
		mycities = GNcities(myPointsUTM,max=5)
		myplot(myPointsUTM)
		
		# and the old school way
		thebbox = bbox(myraster)
		north = thebbox[2,2]
		south = thebbox[2,1]
		east = thebbox[1,2]
		west = thebbox[1,1]
		upperleft = c(north, west)
		lowerright = c(south, east)
		
		mytiles = openmap(upperleft, lowerright, minNumTiles=2)
		mycities = GNcities(north,east,south,west,max=5)
		mycities = spTransform(mycities, CRS(proj4string(myPointsMercator)))
		myplot(myPointsMercator)
		
	}
}