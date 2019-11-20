
earthFile = Pmisc::downloadIfOld(
	'https://www.naturalearthdata.com/http//www.naturalearthdata.com/download/110m/cultural/ne_110m_admin_0_countries_lakes.zip',
	path = tempdir())
earth = raster::shapefile(grep("shp$", earthFile,
	value=TRUE))
earth@data = earth@data[,'NAME', drop=FALSE]


#earth@polygons[[160]]@Polygons[[8]]@coords[,1] = 
#	pmin(pmax(earth@polygons[[160]]@Polygons[[8]]@coords[,1],
#		-179., 179.9)
earth@polygons[[160]]@Polygons[[8]]@coords[,2] = 
	pmax(earth@polygons[[160]]@Polygons[[8]]@coords[,2],
		-87)

worldMap = spTransform(earth, crsMerc)


	myCrsO = moll(worldMap[
		grep("Japan", worldMap$NAME),], 
		angle=25)
 xTcrop = wrapPoly(x=worldMap[15, ], crs=myCrsO)


worldMap = rgeos::gSimplify(
 	worldMap, tol=1000)

 xTcrop = wrapPoly(x=worldMap, crs=myCrsO)


worldMap = SpatialPolygonsDataFrame(
	worldMap, data=earth@data)

save(list='worldMap', 
	file=file.path('/home/patrick/workspace',
		'diseasemapping/pkg/mapmisc/data/worldMap.RData'),
		compress='xz')
