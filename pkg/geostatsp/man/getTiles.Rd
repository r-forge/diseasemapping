\name{getTiles}
\alias{getTiles}
\alias{getTiles.default}
\alias{getTiles.matrix}
\alias{getTiles.Extent}

\title{Get tiles from a server}
\description{
	This function is a wrapper for the 
	\code{\link[webmaps]{getTiles}} 
function in the \link[webmaps]{webmaps} package, where
the limits can be passed as a bounding box or an Extent.
}
\usage{
	getTiles(xlim,  ...)
	\method{getTiles}{matrix}(xlim,...)
	\method{getTiles}{Extent}(xlim,...)
	}

\arguments{
	\item{xlim}{
	A bounding box from a  
	SpatialPoints or SpatialPolygons object or an 
	Extent from a Raster.}
\item{...}{ additional arguments, 
see \code{\link[webmaps]{getTiles}} in the webmaps package.}
}



\examples{

myraster = raster(nrows=11,ncols=11,xmn=0,xmx=10,ymn=0,ymx=10)

\dontrun{
mytiles = getTiles(myraster@extent, zoom=6,path = "http://tile.openstreetmap.org/")

image(mytiles)

mypoints = SpatialPoints(10*cbind(runif(10), runif(10)))

mytiles=getTiles(mypoints@bbox, zoom=6,path = "http://tile.openstreetmap.org/")

image(mytiles)
}


}