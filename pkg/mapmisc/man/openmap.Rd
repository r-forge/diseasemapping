\name{openmap}
\alias{openmap}
\alias{osmTiles}
\alias{openmapAttribution}
\title{Download map tiles}
\description{
Downloads map tiles from Openstreetmap.org and other servers.
}
\usage{
	openmap(x, zoom, 
	path="http://tile.openstreetmap.org/",
	maxTiles = 9,
	crs=projection(x),   buffer=0, 
	fact=1, verbose=FALSE)
	osmTiles(name)
	openmapAttribution(name, type=c('text','latex','markdown','html'), short=FALSE)
}

\arguments{
	\item{x}{
	An \code{\link[raster]{extent}} or any spatial object (raster, Spatial*)
	from which an extent can be obtained.  }
\item{zoom	}{
the zoom level, when missing it will be determined by maxTiles.}
\item{path}{	
the tile server from which to get the map, see \url{http://wiki.openstreetmap.org/wiki/Tiles#Servers}.}
\item{maxTiles}{	
If zoom is missing, zoom will be chosen such that the 
number of map tiles is less than or equl to this number.}
\item{crs}{Projection for the output, defaulting to the same projection as
\code{x}.  If \code{x} has no projection, 
for instance when \code{x} is a matrix or extent, 
\code{crs} is also used as the projection of \code{x}.  If \code{crs} is missing
and \code{x} has no crs, long-lat is used.}
\item{buffer}{Extend the extent for which the map is requested, in units
of \code{x}.  Can be negative, or a vector of length 2 for different x and y 
extensions}
\item{fact}{Passed to \code{\link[raster]{disaggregate}} before 
	reprojecting if \code{fact>1}, 
	helps to produce a clearer image.}
  \item{verbose}{
    If TRUE, give information on where tiles are coming from, cache
    hits, etc
  }
  \item{name}{name of a tile path, if missing a vector of all available tile paths
  is returned.  \code{name} can be any of the names of the vector returned when 
  \code{name} is unspecified.}
  \item{type}{format for the attribution}
  \item{short}{short or long attribution}
}

\details{
These functions download, display, and manipulate map tiles stored in a
standard way either on a web server or a local folder.

Map tiles are a set of PNG images that span the world at a set of zoom
levels. Zoom level 1 has four 256x256 pixel tiles
in a 2x2 pattern over the whole world. In general, zoom level n has 2^n
by 2^n tiles. Zoom levels go up to about 17 or 18 depending on the tile
server.


Be sure to attribute any maps you publish, the osmAttribution function will assist.
}


\value{
	\code{openmap} returns a RasterBrick 
	\code{\link[raster]{brick}}, with 'red', 'green' and 'blue' layers. 
}	

\seealso{\link[OpenStreetMap]{openmap} }

\examples{


myraster = raster(matrix(0,10,10),xmn=8,xmx=18,ymn=0,ymx=10, 
  crs="+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0")
values(myraster) = seq(0,1,len=ncell(myraster))

myPoints = SpatialPoints(myraster, proj4string=CRS(proj4string(myraster)))[
	seq(1,ncell(myraster),len=5)]

names(osmTiles())

\dontrun{

	mytiles = openmap(myraster)

	map.new(myraster)
	plot(mytiles, add=TRUE)
	points(myPoints,col='red')



	mytiles = openmap(myPoints,
		path=osmTiles("mapquest-sat"), verbose=TRUE)
	map.new(myPoints)
	plotRGB(mytiles, add=TRUE)
	points(myPoints,col='red')
	openmapAttribution(mytiles)
}

openmapAttribution("osm", type='markdown')

}