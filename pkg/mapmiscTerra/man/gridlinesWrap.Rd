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
A character string representing a CRS
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

\value{A list with elements \code{lines}, containg the graticule lines, and \code{points} containing the locations and labels for longitude and latitude values.}

\author{
Patrick Brown
}

\seealso{\code{\link[terra]{graticule}}}

\examples{
data('worldMap')
worldMap = terra::unwrap(worldMap)
crsMoll = moll(-100)
worldMapT = wrapPoly(worldMap, crsMoll, buffer.width=200*1000)

plot(attributes(crsMoll)$ellipse)
plot(worldMapT, add=TRUE)
gridlinesWrap(crsMoll, lty=3, col='red', cex=0.6)


}

