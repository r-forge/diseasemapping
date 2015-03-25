insetMap = function(crs, 
    pos="bottomright",map="osm",zoom=0, 
    width=max(c(0.2, 1-par('plt')[2])), 
col="#FF000090", borderSmall=NA, borderBig=NULL,
		cropInset = extent(-180,xmax=180, ymin=-47, ymax=71),
		outer=TRUE) {

fromEdge = matrix(par("plt"), 2, 2, 
		dimnames=list(c("min","max"), c("x","y")))
extentUsr = matrix(par("usr"),2,2, dimnames=dimnames(fromEdge))
dimUsr = abs(apply(extentUsr, 2, diff))
fracUsr = abs(apply(fromEdge, 2, diff))
dimFull = dimUsr/fracUsr

extentFull = extentUsr
extentFull[1,] = extentFull[1,] - dimFull*fromEdge[1,]
extentFull[2,] = extentFull[2,] + dimFull*(1-fromEdge[2,])



if(outer) {
	extentBig = extentFull
}

extentSmall = try(extent(crs), silent=TRUE)

if(class(extentSmall)=="try-error")
	extentSmall = extentBig
	
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

crsCrop = try(CRS(proj4string(cropInset)),silent=TRUE)
if(class(crsCrop)=="try-error")
	crsCrop = crsLL
tocrop = t(bbox(extent(cropInset)))
tocrop = SpatialPoints(tocrop,
		proj4string=crsCrop)

if(is.character(map)) {
  map = openmap(xsp, path=map, zoom=zoom,crs=crsLL)
}
# make sure map is a raster
if(!length(grep("^Raster", class(map)))) {
  warning('map is not a Raster')
}

if(requireNamespace('rgdal', quietly=TRUE)) {
	tocrop = spTransform(tocrop, CRSobj=CRS(proj4string(map)))
	map = crop(map, extent(tocrop))
}

oldinsetbox = t(bbox(map))
oldrange = apply(oldinsetbox, 2, diff)
oldYoverX = oldrange[2]/oldrange[1]

newxrange = diff(par("usr")[1:2])*width

plotFracYcoords = exp(diff(log(apply(matrix(par("usr"),2),2,diff))))

plotFracYinches= exp(diff(log(par('pin'))))

if(length(grep("longlat", projection(map)))) {
	cellRatio = res(area(map))
	cellRatio = cellRatio[1]/cellRatio[2]
} else {
	cellRatio = 1
}
newyrange = newxrange * cellRatio* oldYoverX * plotFracYcoords / plotFracYinches 


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
xsp = coordinates(xsp)

xsp = (xsp - bbOrig[rep(1,N),]) * matrix(scale, N, 2, byrow=TRUE) + 
		bbSmall[rep(1,N),]

if(outer) {
	oldxpd = par("xpd")
	par(xpd=TRUE)
}
if(nlayers(map)>=3) {
 plotRGB(map[[1:3]], add=TRUE)
} else {
	plot(map, add=TRUE)
}
bigpoly = t(bbox(map))
bigpoly = cbind(bigpoly[c(1,2,2,1),1], bigpoly[c(1,1,2,2),2])
polygon(bigpoly,border=borderBig)
 
delta=0.3
theX = anX = c(-delta + delta*1i, -delta + 1i, delta+1i, delta + delta*1i)
for(D in 1:3)
	theX = c(theX, anX*exp(-D*2*pi*1i/4))
theX = theX*exp(-2*pi*1i/8)

if( (diff(range(xsp[,1])))  < (width*dimFull[1]/20) ) {	
	
	polygon((1.5*width*dimFull[1]/20) * theX +
					mean(xsp[,1])+1i*mean(xsp[,2]), col=col,border=NA)
} else {
	polygon(xsp, col=col,border=borderSmall)
}

if(outer) {
	par(xpd=oldxpd)
}	

return(invisible(mapOrig))
}