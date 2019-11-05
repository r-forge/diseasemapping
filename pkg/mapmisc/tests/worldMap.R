#+ setup
library('mapmisc')

if(!exists('mapmiscCachePath'))
  mapmiscCachePath = system.file('extdata', package='mapmisc')

if(!exists('mapmiscCacheReadOnly'))
  mapmiscCacheReadOnly = TRUE

mapmiscCachePath

options(
  mapmiscCachePath = mapmiscCachePath,
  mapmiscCacheReadOnly = mapmiscCacheReadOnly,
  mapmiscVerbose=TRUE)

getOption("mapmiscCachePath")
getOption("mapmiscCacheReadOnly")
#'


#+ themaps, fig.cap='some maps', fig.subcap = c('projectable region', 'projected, with bg','projected, with world countries','projectable madagascar','madagascar')
Spackages = c('rgdal', 'rgeos', 'geosphere')

if(all(unlist(mapply(requireNamespace, package=Spackages, MoreArgs=list(quietly=TRUE))))){
#	library('rgdal')

	data("worldMap")

	world = spTransform(worldMap, crsLL)
	country='Japan'
	Dcountry = grep(country, world$NAME)
	x=world[Dcountry,]

	myCrsO = moll(x, angle=25)

 plot(world)
 plot(attributes(myCrsO)$regionLL, 
 	border='red', col='#0000FF10', add=TRUE)
 
 myMap = openmap(myCrsO, zoom=1, fact=2)
 map.new(myMap)
 plot(myMap, add=TRUE)
 gridlinesWrap(myCrsO, lty=2, col='orange')

 xTcrop = wrapPoly(x=world, crs=myCrsO)

	DcountryT  = grep(country, xTcrop$NAME)
	
	map.new(xTcrop, buffer=1000*1000)
 plot(attributes(myCrsO)$ellipse, add=TRUE, col='lightBlue', border='blue')
	plot(xTcrop,add=TRUE, col='grey')
	plot(xTcrop[DcountryT,], col='red', lty=0, add=TRUE)
	
	gridlinesWrap(myCrsO, lty=2, col='orange')
	

 country='Madagascar'
	Dcountry  = grep(country, world$NAME)
	x=world[Dcountry,]
	
	myCrsMoll = moll(x, angle=100, flip= 'nwu' )
	
	plot(attributes(myCrsMoll)$ellipse)
	
	plot(world)
	plot(attributes(myCrsMoll)$crop, border='red', col='#0000FF10', add=TRUE)
	
	xTcrop = wrapPoly(x=world, crs=myCrsMoll)
	DcountryT  = grep(country, xTcrop$NAME)
	
	map.new(xTcrop)
	plot(attributes(myCrsMoll)$ellipse, add=TRUE, col='lightBlue', border='blue')
	plot(xTcrop,add=TRUE, col='grey')
	plot(xTcrop[DcountryT,], col='green', add=TRUE)
	
	gridlinesWrap(crs=xTcrop, lty=2, col='red')

	
	
	country='Iceland'
	Dcountry  = grep(country, world$NAME)
	x=world[Dcountry,]
	
	myCrsMoll = moll(x,  angle=0, flip= 'nwu' )
	
	plot(attributes(myCrsMoll)$ellipse)
	
	plot(world)
	plot(attributes(myCrsMoll)$crop, border='red', col='#0000FF10', add=TRUE)
	
	xTcrop = wrapPoly(x=world, crs=myCrsMoll)
	DcountryT  = grep(country, xTcrop$NAME)
	
	map.new(xTcrop)
	plot(attributes(myCrsMoll)$ellipse, add=TRUE, col='lightBlue', border='blue')
	plot(xTcrop,add=TRUE, col='grey')
	plot(xTcrop[DcountryT,], col='green', add=TRUE)
	
	gridlinesWrap(crs=xTcrop, lty=2, col='red')
	
	
}
#'