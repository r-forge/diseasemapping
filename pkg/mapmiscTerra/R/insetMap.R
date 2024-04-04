insetMap = function(crs, 
    pos="bottomright",
    map="osm",
    zoom=0, 
    width=max(c(0.2, 1-par('plt')[2])), 
col="#FF000090", borderMap=NULL, 
		cropInset = terra::ext(-180,180, -47, 71),
		outer=TRUE, inset = c(0.1, 0.1),
		...) {

  
if(outer) {
	oldxpd = par("xpd")
	usr = par('usr')
	graphics::clip(usr[1], usr[2], usr[3], usr[4])
	oldpar = par(no.readonly = TRUE) 
	par(xpd=NA)
	on.exit(par(oldpar))    
}
fromEdge = matrix(pmax(0.01, par("plt")), 2, 2, 
		dimnames=list(c("min","max"), c("x","y")))
extentUsr = matrix(par("usr"),2,2, dimnames=dimnames(fromEdge))
dimUsr = abs(apply(extentUsr, 2, diff))
fracUsr = abs(apply(fromEdge, 2, diff))
dimFull = dimUsr/fracUsr

extentFull = extentUsr
extentFull['max',] = extentFull['min',] +
    apply(extentUsr,2,diff) / (fromEdge['max',])
#extentFull['min',] = extentFull['min',] -
#    apply(extentUsr,2,diff) *(  fromEdge['min',])/ (1 -  fromEdge['min',])

extentFull = terra::ext(
		extentFull[1,1],
		extentFull[2,1],
		extentFull[1,2],
		extentFull[2,2]
		)

extentUsr = terra::ext(
		extentUsr[1,1],
		extentUsr[2,1],
		extentUsr[1,2],
		extentUsr[2,2]
)
		

if(outer) {
	extentBig = extentFull
} else {
	extentBig = extentUsr
}

extentSmall = extentUsr
		

bboxSmall = matrix(as.vector(extentSmall), 2)

xseq = seq(extentSmall$xmin, extentSmall$xmax,len=20)
yseq = seq(extentSmall$ymin, extentSmall$ymax,len=20)

polySmall = cbind( 
		c(xseq, rep(extentSmall$xmax, length(yseq)), 
				rev(xseq), rep(extentSmall$xmin, length(yseq))), 
		c(rep(extentSmall$ymin, length(xseq)), yseq,
				rep(extentSmall$ymax, length(xseq)), rev(yseq)
		)
)

crs = terra::crs(crs)
xsp = vect(polySmall, 	crs = crs)


# if cropInset is numeric
  # use it to extend the extent of the plot region
  # and crop the inset map
	
if(is.numeric(cropInset)) {
  cropInset = rast(nrows=40, ncols=40, extent=terra::extend(extentSmall, cropInset), crs=crs)
}

	mapExtent = xsp

if(all(class(cropInset)=='SpatExtent')){
	cropInset = rast(cropInset, crs=crsLL)
} 
# map = 'osm'
if(is.character(map)) {
	cropInsetVec = terra::vect(terra::xyFromCell(cropInset, 1:terra::ncell(cropInset)), crs=crs(cropInset))
	forMap = terra::union(project(mapExtent, crsMerc), project(cropInsetVec, crsMerc))
  map = openmap(x=forMap, path=map, zoom=zoom,crs=crsMerc, verbose=TRUE)
  mapOrig = terra::deepcopy(map)
} else {
	map = terra::deepcopy(map)
	mapOrig = terra::deepcopy(map)
}


# make sure map is a raster

if(!length(grep("^SpatRaster", class(map)))) {
  warning('map is not a Raster')
}

if(!is.null(cropInset)) {
tocrop = terra::rast(
  terra::union(
  	terra::ext(cropInset), 
    terra::ext(terra::project(terra::ext(mapExtent), terra::crs(mapExtent), terra::crs(cropInset)))
  ), crs=terra::crs(cropInset))
tocrop = terra::project(tocrop,   terra::crs(map))
# if the extents are overlapping, crop
map = terra::crop(map, ext(tocrop))
mapOrig = terra::deepcopy(map)
}



usr = par('usr')

plotCellRatio = diff(usr)[-2]/par('pin')
mapSize = diff(as.vector(terra::ext(map)))[-2]

mapXwidth = diff(usr[1:2])*width
mapYwidth = mapXwidth * (plotCellRatio[2]/plotCellRatio[1]) * (mapSize[2]/mapSize[1])
mapYwidthUsr = mapXwidth  * (mapSize[2]/mapSize[1])


if(is.character(pos)) {

	xpoints = matrix(as.vector(extentBig), ncol=2)

	x = apply(xpoints, 2, mean) - 0.5*c(mapXwidth, mapYwidth)


	theRange = diff(usr)[-2]
	inset = rep(inset, 2)

if(length(grep("^top",pos)))
	x[2] = usr[4]-theRange[2]*inset[2] - mapYwidth
if(length(grep("^bottom",pos)))
	x[2] = usr[3]+theRange[2]*inset[2]
if(length(grep("right$",pos)))
	x[1] = usr[2]-theRange[1]*inset[1] - mapXwidth
if(length(grep("left$",pos)))
	x[1] = usr[1]+theRange[1]*inset[1]
} else {
	x=pos
}




terra::ext(map)= terra::ext(c(x[1], x[1]+mapXwidth, x[2], x[2]+mapYwidth))
terra::crs(map) = crs

bbOrig = matrix(as.vector(terra::ext(mapOrig)), 2)
bbSmall = matrix(as.vector(terra::ext(map)), 2)



xsp = terra::project(xsp, terra::crs(mapOrig))
xsp = terra::crop(xsp, mapOrig)

scale =  apply(bbSmall, 2, diff)/ apply(bbOrig, 2, diff)

N = length(xsp)


if(N) {
xsp = vect(
	(terra::crds(xsp) - bbOrig[rep(1,N),]) * matrix(scale, N, 2, byrow=TRUE) + bbSmall[rep(1,N),], 
	crs = crs(mapOrig))
}

toScale = list(shift1 = bbOrig[1,], scale=scale, shift2 = bbSmall[1,])		

xsp = terra::crds(xsp)

if(terra::has.RGB(map)) {
	theCol = do.call(grDevices::rgb, 
		c(as.list(as.data.frame(values(map)[,c('red','green','blue')])), 
			list(maxColorValue=255))
	) 
	theIndex = matrix(1:length(theCol), ncol(map), nrow(map))[,nrow(map):1]
} else {
	if(all(c('red','green','blue') %in% names(terra::coltab(map)[[1]]))) {
		theCol = do.call(grDevices::rgb, 
			c(as.list(terra::coltab(map)[[1]][,c('red','green','blue')]), list(maxColorValue=255)))	
		theIndex = 	matrix(terra::values(map), ncol(map), nrow(map) )[, nrow(map):1]
	} else {
		theCol = NA
		theIndex = matrix(0, ncol(map), nrow(map))
	}

}
do.call(graphics::clip, as.list(par('usr')))
graphics::image(
	terra::xFromCol(map), 
	terra::yFromRow(map, nrow(map):1), 
	theIndex,
	useRaster=FALSE, add=TRUE, col=theCol)
# for some reason, the code below crops the map
#plot(map, add=TRUE)

# border around the map
bigpoly = matrix(as.vector(terra::ext(map)), ncol=2)
bigpoly = cbind(bigpoly[c(1,2,2,1),1], bigpoly[c(1,1,2,2),2])

graphics::polygon(bigpoly,border=borderMap)
 
delta=0.3
theX = anX = c(-delta + delta*1i, -delta + 1i, delta+1i, delta + delta*1i)
for(D in 1:3)
	theX = c(theX, anX*exp(-D*2*pi*1i/4))
theX = theX*exp(-2*pi*1i/8)

if( (diff(range(xsp[,1])))  < (width*dimFull[1]/2000) ) {	
	graphics::polygon((1.5*width*dimFull[1]/20) * theX +
					mean(xsp[,1])+1i*mean(xsp[,2]), col=col, ...)
} else {
	graphics::polygon(xsp, col=col, border=substr(col, 1, 7), ...)
}



attributes(xsp)$toScale = toScale

return(invisible(xsp))
}