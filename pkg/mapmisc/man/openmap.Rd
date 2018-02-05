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
	  crs=projection(x),   
	  buffer=0, fact=1, 
	 verbose=getOption('mapmiscVerbose'),
	  cachePath=getOption('mapmiscCachePath')
)

	osmTiles(name, xyz, suffix)

	openmapAttribution(name, 
	  type=c('text','latex','markdown','html', 'auto'), 
	  short=FALSE)
}

\arguments{
	\item{x}{
	An \code{\link[raster]{extent}} or any spatial object (raster, Spatial*)
	from which an extent can be obtained.  }
\item{zoom	}{
the zoom level, when missing it will be determined by maxTiles.}
\item{path}{	
Source of map tiles, see \url{http://diseasemapping.r-forge.r-project.org/mapLayers.html}.}
\item{maxTiles}{	
If zoom is missing, zoom will be chosen such that the 
number of map tiles is less than or equl to this number.}
\item{crs}{Projection for the output, defaulting to the same projection as
\code{x}.  If \code{x} has no projection, 
for instance when \code{x} is a matrix or extent, 
\code{crs} is also used as the projection of \code{x}.  If \code{crs} is \code{NA} or missing and \code{x} has no crs, long-lat is used.  }
\item{buffer}{Extend the extent for which the map is requested, in units
of \code{x}.  Can be negative, or a vector of length 2 for different x and y 
extensions}
\item{fact}{Passed to \code{\link[raster]{disaggregate}} before 
	reprojecting if \code{fact>1}, 
	helps to produce a clearer image.}
  \item{verbose}{
    Print information about map images being downloaded, defaults to \code{FALSE}.
  }
  \item{cachePath}{Location to store downloaded map images, defaults to \code{tempdir()}}
  \item{name}{name of a tile path, if missing a vector of all available tile paths
  is returned.  \code{name} can be any of the names of the vector returned when 
  \code{name} is unspecified.}
  \item{type}{format for the attribution}
  \item{short}{short or long attribution}
  \item{xyz}{format of xyz coordinates in URL's}
  \item{suffix}{string to append to URL's, i.e. \code{'.png'}}
}

\details{
These functions download, display, and manipulate map tiles stored in a
standard way either on a web server or a local folder.

Map tiles are a set of PNG images that span the world at a set of zoom
levels. Zoom level 1 has four 256x256 pixel tiles
in a 2x2 pattern over the whole world. In general, zoom level n has \eqn{2^n}
by \eqn{2^n} tiles. Zoom levels go up to about 17 or 18 depending on the tile
server.


See \url{http://diseasemapping.r-forge.r-project.org/mapLayers.html} for a full set of 
map tiles.

Be sure to attribute any maps you publish, the \code{osmAttribution} function will assist.  If \code{type = 'auto'} then markdown format will be used unless a variable \code{mdToTex} is defined and equal to \code{TRUE}.

}


\value{
	\code{openmap} returns a \code{Raster} with indexed colours,
	 or a \code{RasterBrick} 
	\code{\link[raster]{brick}}, with 'red', 'green' and 'blue' layers.
	
	\code{openmapAttribution} returns a character string.
	 
}	

\seealso{\link[OpenStreetMap]{openmap}}

\examples{

data("netherlands")

plot(nldTiles)

openmapAttribution(nldTiles, short=TRUE, type='markdown')

openmapAttribution("stamen-toner", type='text')


myraster = raster(matrix(0,10,10),xmn=8,xmx=18,ymn=0,ymx=10, 
  crs="+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0")

values(myraster) = seq(0,1,len=ncell(myraster))

myPoints = SpatialPoints(myraster, proj4string=CRS(proj4string(myraster)))[
	seq(1,ncell(myraster),len=5)]

names(osmTiles())


\dontrun{
	mytiles = openmap(myraster, zoom=6)

	map.new(myraster)
	plot(mytiles, add=TRUE)
	points(myPoints,col='red')

	mytiles = openmap(myPoints,
		path='waze', verbose=TRUE)
	map.new(myPoints)
	plotRGB(mytiles, add=TRUE)
	points(myPoints, col='red')
	openmapAttribution(mytiles)
}

}