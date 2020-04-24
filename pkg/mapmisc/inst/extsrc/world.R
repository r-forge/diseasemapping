library('mapmisc')
earthFile = Pmisc::downloadIfOld(
	'https://www.naturalearthdata.com/http//www.naturalearthdata.com/download/110m/cultural/ne_110m_admin_0_countries_lakes.zip',
	path = tempdir())
earth = raster::shapefile(grep("shp$", earthFile,
	value=TRUE))
earth@data = earth@data[,'NAME', drop=FALSE]
earth$NAME = gsub('^C..te ', 'Cote ', earth$NAME)

earth2 = rgeos::gSimplify(earth, tol=0.001, TRUE)

#earth@polygons[[160]]@Polygons[[8]]@coords[,1] = 
#	pmin(pmax(earth@polygons[[160]]@Polygons[[8]]@coords[,1],
#		-179., 179.9)
earth2@polygons[[160]]@Polygons[[8]]@coords[,2] = 
	pmax(earth2@polygons[[160]]@Polygons[[8]]@coords[,2],
		-87)

stuff = list()
for(D in 1:nrow(earth@data)) {
  cat(D, ' ')
  stuff[[D]] = try(spTransform(earth2[D,], crsMerc))
}
unlist(lapply(stuff, class))

worldMap = spTransform(earth2[1:nrow(earth@data)], crsMerc)
worldMap = SpatialPolygonsDataFrame(worldMap, earth@data)

	myCrsO = moll(worldMap[
		grep("Japan", worldMap$NAME),])
	
	plot(attributes(myCrsO)$regionLL)
bbox(attributes(myCrsO)$regionLL)
	
 xTcrop = try(wrapPoly(x=worldMap[15, ], crs=myCrsO))


worldMap = rgeos::gSimplify(
 	worldMap, tol=1000)

 xTcrop = wrapPoly(x=worldMap, crs=myCrsO)


worldMap = SpatialPolygonsDataFrame(
	worldMap, data=earth@data)

worldMap$NAME = gsub("C.te d.Ivoire", "Cote d Ivoire", worldMap$NAME)


save(list='worldMap', 
	file=file.path('~',
#	  'workspace',
    'research',
		'diseasemapping/pkg/mapmisc/data/worldMap.RData'),
		compress='xz')
