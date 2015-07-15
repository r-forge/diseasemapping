\name{map.new}
\alias{map.new}


\title{
Start a new map
}
\description{
Prepare a plotting window suitable for a map
}
\usage{
map.new(x,legendRight=FALSE, buffer=0, mar=c(0,0,0,0), ...)
}

\arguments{
\item{x}{
A spatial object from which an extent can be extracted.
}
\item{legendRight}{Leave room to the right for 
the legend produced by plotting a Raster object}
\item{buffer}{passed to \code{\link[raster]{extend}} to increase the plotting area}
\item{mar}{see \code{\link[graphics]{par}}}
\item{\dots}{
Additional arguments passed to \code{plot} 
}
}
\details{
\code{map.new} initiates a plot intended to contain a map covering the extent of \code{x}, 
with no margins.
}
\author{
Patrick Brown
}

\seealso{\code{\link[raster]{scalebar}}, \code{\link[sp]{spplot}}}

\examples{

Npoints = 20
myPoints = SpatialPointsDataFrame(
	cbind(runif(Npoints), 51+runif(Npoints)), 
	data=data.frame(y1=c(NA, rnorm(Npoints-1)), 
	y2=c(sample(0:5, Npoints-1,replace=TRUE), NA)),
	proj4string=CRS(
	"+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0"
	)
)

map.new(myPoints, legendRight=TRUE, mar=c(3,3,0,0), buffer=0.2)
points(myPoints)
legendBreaks('right', list(breaks=1:3, col=1:2))

}

