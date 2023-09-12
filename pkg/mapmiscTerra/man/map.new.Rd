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
\item{buffer}{passed to \code{\link[terra]{extend}} to increase the plotting area}
\item{mar}{see \code{\link[graphics]{par}}}
\item{\dots}{
Additional arguments passed to \code{plot} 
}
}
\details{
\code{map.new} initiates a plot intended to contain a map covering the extent of \code{x}, 
with no margins.
}

\value{A list of the graphical parameters prior to calling \code{map.new}}
\author{
Patrick Brown
}

%\seealso{}

\examples{

nldTiles = terra::unwrap(nldTiles)
nldCities = terra::unwrap(nldCities)


oldpar = map.new(nldCities)
plot(nldTiles, add=TRUE)
points(nldCities)
par(oldpar)

}

