#Produced by polyHédronisme http://levskaya.github.com/polyhedronisme
# group dccdkdkdkD

x = read.table(file.path(Sys.getenv("HOME"), 'Downloads', 
	'polyhedronisme-A100dccdktkD.obj'),
	skip = 3, header=FALSE, stringsAsFactors=FALSE,
	nrow=4322)


x2 = as.matrix(x[,-1])
rgl::plot3d(x2)

polyhedron1 = cbind(
	atan(x2[,2]/x2[,1]),
	acos(x2[,3]/sqrt(apply(x2^2, 1, sum)))
	)

plot(polyhedron1)

polyhedron2 = cbind(
	polyhedron1[,1]*360/pi,
	polyhedron1[,2]*90*2/pi-90
	)
plot(polyhedron2)

polyhedron = sp::SpatialPoints(polyhedron2, proj4string=mapmisc::crsLL)

myMap = mapmisc::openmap(polyhedron)

raster::plot(myMap)
sp::plot(polyhedron, add=TRUE)


dput(polyhedron)




  llPoints = polyhedron@coords
  
  Sprob = seq(0,1,len=251)
  latSeq = sort(unique(c(
        seq(-90,90,len=201),
        quantile(llPoints[,2], prob = Sprob))))
  lonSeq = seq(-180,180,len=51)
  
  llBorder = sp::SpatialPoints(cbind(
    lon=c(
      lonSeq, rep(-180, length(latSeq)),
      rep(180, length(latSeq))
    ),
    lat=c(
      rep(-90,length(lonSeq)), latSeq, latSeq
    )), proj4string=mapmisc::crsLL)



  raster::plot(myMap)
	sp::plot(llBorder, add=TRUE)

dput(llBorder)

saveRDS(llBorder, ascii=TRUE, compress=FALSE, file='stuff')
saveRDS(polyhedron, ascii=TRUE, compress=TRUE, file='stuff')



myMap2 = projectRaster(myMap, crs=milliondeaths::indiaCrs, method='ngb')
myMap2@legend@colortable = myMap@legend@colortable

plot(myMap2)
