insetMap = function(crs, 
    pos="bottomright",map="osm",zoom=1, 
    width=max(c(0.2, 1-par('plt')[2])), 
col="#FF000090", borderMap=NULL, 
		cropInset = terra::ext(-180,180, -47, 71),
		outer=TRUE, ...) {

  
  if(outer) {
	oldxpd = par("xpd")
	par(xpd=TRUE)
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
  cropInset = rast(terra::extend(extentSmall, cropInset), crs=crs)
}

if(all(class(cropInset)=='SpatExtent')){
	cropInset = rast(cropInset, crs=crsLL)
	mapExtent = xsp
} else {
	if(is.null(cropInset)) {
		mapExtent = xsp
	}	 else {
		mapExtent = cropInset
	}
}

if(is.character(map)) {
  map = openmap(x=mapExtent, path=map, zoom=zoom,crs=crsMerc, verbose=TRUE)
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
tocrop = terra::project(terra::ext(tocrop), terra::crs(tocrop), terra::crs(map))
# if the extents are overlapping, crop
map = terra::crop(map, tocrop)
mapOrig = terra::deepcopy(map)
}





oldrange = diff(as.vector(extentFull))[-2]
oldYoverX = oldrange[2]/oldrange[1]

newxrange = (terra::xmax(extentFull)-terra::xmin(extentFull))*width
plotFracYcoords = oldrange['y']/oldrange['x']
plotFracYinches= par('pin')[2]/par('pin')[1]
    

cellRatio = apply(matrix(par('usr'), 2),2,diff)/par('pin')

insetMapRatio = diff(as.vector(terra::ext(map)))[-2]
insetMapRatio = insetMapRatio[2]/insetMapRatio[1]
  
newyrange = newxrange * insetMapRatio# * plotFracYcoords / plotFracYinches 


if(is.character(pos)) {

	xpoints = matrix(as.vector(extentBig), ncol=2)

	x = apply(xpoints, 2, mean) - 0.5*c(newxrange, newyrange)

if(length(grep("^top",pos)))
	x[2] = xpoints[2,2]-newyrange
if(length(grep("^bottom",pos)))
	x[2] = xpoints[1,2]
if(length(grep("right$",pos)))
	x[1] = xpoints[2,1]-newxrange
if(length(grep("left$",pos)))
	x[1] = xpoints[1,1]
} else {
	x=pos
}



terra::ext(map)= terra::ext(c(x[1], x[1]+newxrange, x[2],
				x[2]+newyrange*(cellRatio[2]/cellRatio[1])))
terra::crs(map) = crs
bbOrig = matrix(as.vector(terra::ext(mapOrig)), 2)
bbSmall = matrix(as.vector(terra::ext(map)), 2)


xsp = terra::project(xsp, terra::crs(mapOrig))

scale =  apply(bbSmall, 2, diff)/ apply(bbOrig, 2, diff)

N = length(xsp)



xsp = vect(
	(terra::crds(xsp) - bbOrig[rep(1,N),]) * matrix(scale, N, 2, byrow=TRUE) + bbSmall[rep(1,N),], 
	crs = crs(mapOrig))

toScale = list(shift1 = bbOrig[1,], scale=scale, shift2 = bbSmall[1,])		

xsp = terra::crop(xsp, mapOrig)
xsp = terra::crds(xsp)



plot(map, add=TRUE)

# border around the map
bigpoly = matrix(as.vector(terra::ext(map)), ncol=2)
bigpoly = cbind(bigpoly[c(1,2,2,1),1], bigpoly[c(1,1,2,2),2])

graphics::polygon(bigpoly,border=borderMap)
 
delta=0.3
theX = anX = c(-delta + delta*1i, -delta + 1i, delta+1i, delta + delta*1i)
for(D in 1:3)
	theX = c(theX, anX*exp(-D*2*pi*1i/4))
theX = theX*exp(-2*pi*1i/8)

if( (diff(range(xsp[,1])))  < (width*dimFull[1]/20) ) {	
	graphics::polygon((1.5*width*dimFull[1]/20) * theX +
					mean(xsp[,1])+1i*mean(xsp[,2]), col=col, ...)
} else {
	graphics::polygon(xsp, col=col, border=NA, ...)
}

if(outer) {
	par(xpd=oldxpd)
}	

attributes(xsp)$toScale = toScale

return(invisible(xsp))
}