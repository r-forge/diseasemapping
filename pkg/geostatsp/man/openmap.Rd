\name{openmap}
\alias{openmap}
\title{Get map from a server}
\description{
	This function is a wrapper for the 
	\code{\link[OpenStreetMap]{openmap}} 
function in the \link[OpenStreetMap]{openmap} package
}
\usage{
	openmap(upperLeft, lowerRight=NULL, zoom = NULL,
		type = c("osm", "osm-bw", "maptoolkit-topo", "waze", "mapquest", "mapquest-aerial", "bing", "stamen-toner", "stamen-terrain", "stamen-watercolor", "osm-german", "osm-wanderreitkarte", "mapbox", "esri", "esri-topo", "nps", "apple-iphoto", "skobbler", "cloudmade-<id>", "hillshade", "opencyclemap", "osm-transport", "osm-public-transport", "osm-bbike", "osm-bbike-german"),
		minNumTiles = 9L, mergeTiles = TRUE) 
}

\arguments{
	\item{upperLeft}{
	Either a vector of lat and long for the lower left corner (note: c(lat,long), not c(long,lat)!),	
a raster, SPDF, bounding box or an Extent.  Bounding boxes and extents must be in long-lat.}
\item{lowerRight}{the lower right lat and long, ignored unless \code{upperLeft} 
is a vector}
\item{
zoom	}{
the zoom level. If null, it is determined automatically.}
\item{
type}{	
the tile server from which to get the map.}
\item{
minNumTiles}{	
If zoom is null, zoom will be chosen such that the number of map tiles is greater than or equal to this number..}
\item{
mergeTiles}{	
should map tiles be merged into one tile.}
}

\value{
	an object of class \code{OpenStreetMap}, in the same 
	projection as \code{upperLeft} if a projection was available, 
	otherwise in a mercator projection
}	

\seealso{\code{\link[OpenStreetMap]{openmap}}}

\examples{


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

\dontrun{
	par(mar=c(0,0,0,0))
	mytiles = openmap(myraster, minNumTiles=2,  type="osm")
	plot(myraster)
	plot(mytiles, add=TRUE)
	points(myPoints,col='red')

	mytiles = openmap(extent(myraster), minNumTiles=2,  type="osm")
	plot(myPointsMercator)
	plot(mytiles, add=TRUE)
	points(myPointsMercator,col='red')

	mytiles = openmap(bbox(myraster), minNumTiles=2,   type="osm")
	plot(myPointsMercator)
	plot(mytiles, add=TRUE)
	points(myPointsMercator,col='red')

	mytiles = openmap(myPoints, minNumTiles=2,  type="osm")
	plot(myPoints)
	plot(mytiles, add=TRUE)
	points(myPoints,col='red')
	

	mytiles = openmap(myrasterUTM, minNumTiles=2,  type="osm")
	plot(myrasterUTM)
	plot(mytiles, add=TRUE)
	points(myPointsUTM,col='red')


	mytiles = openmap(myPointsUTM, minNumTiles=2,  type="osm")
	plot(myPointsUTM)
	plot(mytiles, add=TRUE)
	points(myPointsUTM,col='red')
	
	ex = extent(myraster)
	mytiles = openmap(c(ex@ymax, ex@xmin), c(ex@ymin, ex@xmax),
		minNumTiles=2, type="osm")
	plot(myPointsMercator)
	plot(mytiles, add=TRUE)
	points(myPointsMercator,col='red')
		
	
}

}