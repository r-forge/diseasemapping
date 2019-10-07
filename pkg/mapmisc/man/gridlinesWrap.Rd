\name{gridlinesWrap}
\alias{gridlinesWrap}

\title{
Adds long-lat grid for projected data
}
\description{
long-lat grid lines are added to a map in the coordinate system specified, 
allowing for map projections wrapped differently from the 180 meridian.
}
\usage{
gridlinesWrap(crs, 
		easts=seq(-180,180,by=60),
		norths=seq(-90,90,by=30),
		ndiscr=40, plotLines=TRUE, 
		plotLabels = TRUE, ...)
}

\arguments{
  \item{crs}{
A \code{CRS} object, proj4 string, or an object from which a projection 
can be extracted with \code{proj4string(crs)}
}
\item{easts}{vector of longitudes}
\item{norths}{vector of latitudes}
\item{ndiscr}{number of intermediate points per line}
\item{plotLines}{add lines to existing plot}
\item{plotLabels}{add labels to existing plot}
\item{\dots}{
Additional arguments passed to \code{lines} or \code{text}, 
for example \code{lty=2} 
}
}

\author{
Patrick Brown
}

\seealso{\code{\link[sp]{gridlines}}, \code{\link[rgdal]{llgrid}}}

\examples{
\dontrun{
Npoints = 20
myPoints = SpatialPointsDataFrame(
	cbind(
	runif(Npoints, -15000000, 15000000), 
	runif(Npoints, -8000000, 8000000)), 
	data=data.frame(y1=c(NA, rnorm(Npoints-1)), 
	y2=c(sample(0:5, Npoints-1,replace=TRUE), NA)),
	proj4string=moll(c(-100,0))
)


plot(myPoints)
gridlinesWrap(myPoints, lty=3, col='red')
}

}

