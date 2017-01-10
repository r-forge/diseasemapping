insetMap = function(crs, 
    pos="bottomright",map="osm",zoom=0, 
    width=max(c(0.2, 1-par('plt')[2])), 
col="#FF000090", borderMap=NULL, 
		cropInset = extent(-180,xmax=180, ymin=-47, ymax=71),
		outer=TRUE, ...) {

  
  
fromEdge = matrix(par("plt"), 2, 2, 
		dimnames=list(c("min","max"), c("x","y")))
extentUsr = matrix(par("usr"),2,2, dimnames=dimnames(fromEdge))
dimUsr = abs(apply(extentUsr, 2, diff))
fracUsr = abs(apply(fromEdge, 2, diff))
dimFull = dimUsr/fracUsr

extentFull = extentUsr
extentFull['max',] = extentFull['min',] +
    apply(extentUsr,2,diff) / (fromEdge['max',])
extentFull['min',] = extentFull['min',] -
    apply(extentUsr,2,diff) *(  fromEdge['min',])/ (1 -  fromEdge['min',])

extentFull = extent(
		extentFull[1,1],
		extentFull[2,1],
		extentFull[1,2],
		extentFull[2,2]
		)

extentUsr = extent(
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
		
if(is.character(crs))
			crs = CRS(crs)
if(class(crs) != "CRS")
	crs = CRS(proj4string(crs))

bboxSmall = t(bbox(extentSmall))

xseq = seq(bboxSmall[1,1], bboxSmall[2,1],len=20)
yseq = seq(bboxSmall[1,2], bboxSmall[2,2],len=20)

polySmall = cbind( 
		c(xseq, rep(bboxSmall[2,1], length(yseq)), 
				rev(xseq), rep(bboxSmall[1,1], length(yseq))), 
		c(rep(bboxSmall[1,2], length(xseq)), yseq,
				rep(bboxSmall[2,2], length(xseq)), rev(yseq)
		)
)


xsp = SpatialPoints(polySmall, 	proj4string = crs)


# if cropInset is numeric
  # use it to extend the extent of the plot region
  # and crop the inset map
	
if(is.numeric(cropInset)) {
  cropInset = raster(raster::extend(extentSmall, cropInset), crs=crs)
}

if(all(class(cropInset)=='Extent')){
	cropInset = raster(cropInset, crs=crsLL)
	mapExtent = xsp
} else {
	if(is.null(cropInset)) {
		mapExtent = xsp
	}	 else {
		mapExtent = cropInset
	}
}

if(is.character(map)) {
  map = openmap(mapExtent, path=map, zoom=zoom,crs=NA)
}
# make sure map is a raster
if(!length(grep("^Raster", class(map)))) {
  warning('map is not a Raster')
}

if(!is.null(cropInset)) {
tocrop = raster(
  raster::union(extent(cropInset), 
    extent(projectExtent(mapExtent, projection(cropInset)))
  ), crs=projection(cropInset))
tocrop = projectExtent(tocrop, CRS(proj4string(map)))
# if the extents are overlapping, crop
map = crop(map, tocrop)
}






oldrange = apply(bbox(extentFull), 2, diff)
oldYoverX = oldrange[2]/oldrange[1]

newxrange = (xmax(extentFull)-xmin(extentFull))*width
plotFracYcoords = oldrange['y']/oldrange['x']
plotFracYinches= par('pin')[2]/par('pin')[1]
    
#plotFracYcoords = exp(diff(log(apply(matrix(par("usr"),2),2,diff))))

#plotFracYinches= exp(diff(log(par('pin'))))

#if(length(grep("longlat", projection(map)))) {
#	cellRatio = res(area(map))
#	cellRatio = cellRatio[1]/cellRatio[2]
#} else {
	cellRatio = 1
#}

insetMapRatio = abs(apply(bbox(map),1,diff))
insetMapRatio = insetMapRatio[2]/insetMapRatio[1]
  
newyrange = newxrange * cellRatio* insetMapRatio# * plotFracYcoords / plotFracYinches 


if(is.character(pos)) {
		xpoints = t(bbox(extentBig))

x = apply(xpoints, 2, mean) - 0.5*c(newxrange, newyrange)

if(length(grep("^top",pos)))
	x[2] = xpoints[2,2]-newyrange
if(length(grep("^bottom",pos)))
	x[2] = xpoints[1,2]
if(length(grep("right$",pos)))
	x[1] = xpoints[2,1]-newxrange
if(length(grep("left$",pos)))
	x[1] = xpoints[1,1]
} else x=pos


mapOrig = map
extent(map)= extent(c(x[1], x[1]+newxrange, x[2],
				x[2]+newyrange))
proj4string(map) = CRS()
bbOrig = t(bbox(extent(mapOrig)))
bbSmall = t(bbox(extent(map)))


if(requireNamespace('rgdal', quietly=TRUE)) {
  xsp = spTransform(xsp,
		CRSobj=CRS(proj4string(mapOrig)))
}

scale =  apply(bbSmall, 2, diff)/ apply(bbOrig, 2, diff)

N = length(xsp)

xsp@coords = (xsp@coords - bbOrig[rep(1,N),]) * matrix(scale, N, 2, byrow=TRUE) + 
		bbSmall[rep(1,N),]

xsp = raster::crop(xsp, map)
xsp = coordinates(xsp)

if(outer) {
	oldxpd = par("xpd")
	par(xpd=TRUE)
}
if(nlayers(map)>=3) {
 plotRGB(map[[1:3]], add=TRUE)
} else {
	plot(map, add=TRUE)
}

# border around the map
bigpoly = t(bbox(map))
bigpoly = cbind(bigpoly[c(1,2,2,1),1], bigpoly[c(1,1,2,2),2])
polygon(bigpoly,border=borderMap)
 
delta=0.3
theX = anX = c(-delta + delta*1i, -delta + 1i, delta+1i, delta + delta*1i)
for(D in 1:3)
	theX = c(theX, anX*exp(-D*2*pi*1i/4))
theX = theX*exp(-2*pi*1i/8)

if( (diff(range(xsp[,1])))  < (width*dimFull[1]/20) ) {	
	polygon((1.5*width*dimFull[1]/20) * theX +
					mean(xsp[,1])+1i*mean(xsp[,2]), col=col, ...)
} else {
	polygon(xsp, col=col, border=NA, ...)
}

if(outer) {
	par(xpd=oldxpd)
}	

return(invisible(xsp))
}