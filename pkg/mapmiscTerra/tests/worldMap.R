#+ setup
options(warn=1)

library('mapmisc')

#if(!exists('mapmiscCachePath'))
#  mapmiscCachePath = system.file('extdata', package='mapmisc')

#if(!exists('mapmiscCacheReadOnly'))
#  mapmiscCacheReadOnly = TRUE

#mapmiscCachePath

#options(
#  mapmiscCachePath = mapmiscCachePath,
#  mapmiscCacheReadOnly = mapmiscCacheReadOnly,
#  mapmiscVerbose=TRUE)

getOption("mapmiscCachePath")
#'

#+ themaps, fig.cap='some maps', fig.subcap = c('projectable region', 'projected, with bg','projected, with world countries','projectable madagascar','madagascar')


		data("worldMap")
		worldMap = unwrap(worldMap)

		world = project(worldMap, crsLL)
		country='Japan'
		x=world[grep(country, world$NAME)]

		myCrsO = moll(x)
if(!interactive()) pdf("worldMap.pdf")

 plot(world, ylim = 90.5*c(-1,1))
 plot(attributes(myCrsO)$crop, 
 	border='red', col='#0000FF10', add=TRUE)

#x = myCrsO;zoom=1;fact=0.7;path='osm';crs = crs(x);buffer=0;verbose=TRUE;cachePath = tempdir();suffix=NULL  

 myMap = openmap(myCrsO, zoom=1, fact=0.7, verbose=TRUE)

 plot(myMap)
 plot(attributes(myCrsO)$ellipse, add=TRUE, border='black')
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
	
	myCrsMoll = moll(x)
	
	plot(world)
	plot(attributes(myCrsMoll)$crop, border='red', col='#0000FF10', add=TRUE)


	xTcrop = wrapPoly(x=world, crs=myCrsMoll)
	plot(attributes(myCrsMoll)$ellipse, col='lightBlue')
	plot(xTcrop, add=TRUE, col='white')	

	DcountryT  = grep(country, xTcrop$NAME)
	
	map.new(xTcrop)
	plot(attributes(myCrsMoll)$ellipse, add=TRUE, col='lightBlue', border='blue')
	plot(xTcrop,add=TRUE, col='grey')
	plot(xTcrop[DcountryT,], col='green', add=TRUE)
	
	gridlinesWrap(crs=myCrsMoll, lty=2, col='red')


	
	country='Iceland'
	Dcountry  = grep(country, world$NAME)
	x=world[Dcountry,]
	
	myCrsMoll = moll(x,  angle=10)
	xTcrop = wrapPoly(x=world, crs=myCrsMoll)
	
	plot(attributes(myCrsMoll)$ellipse, col='lightBlue')
	plot(xTcrop, add=TRUE, col='white')
	
	plot(world)
	plot(attributes(myCrsMoll)$crop, border='red', col='#0000FF10', add=TRUE)
	
	DcountryT  = grep(country, xTcrop$NAME)
	
	map.new(myCrsMoll, col='lightBlue')
	plot(xTcrop,add=TRUE, col='grey')
	gridlinesWrap(crs=xTcrop, lty=2, col='red')
	plot(xTcrop[DcountryT,], col='green', add=TRUE)


	# check openmap
	myMap = openmap(myCrsMoll)
	plot(myMap)
	
	
	if(!interactive()) dev.off()
	
#'


