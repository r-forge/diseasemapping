\name{scaleBar}
\alias{scaleBar}
\alias{map.new}

\title{
Map with scale bar
}
\description{
Utilities for plotting a map, adding a scale bar and north arrow, 
and adding a legend of colour scales.
}
\usage{
scaleBar(crs, pos = "bottomright", scale.cex = 1, ...)
map.new(x,rasterLegend=FALSE)
}

\arguments{
  \item{crs}{
A \code{CRS} object, proj4 string, or an object from which a projection 
can be extracted with \code{proj4string(crs)}
}
  \item{pos}{
Position, as specified in the \code{legend} function.
}
  \item{scale.cex}{
Scaling factor for the scale bar.
}
\item{x}{
A spatial object from which an extent can be extracted.
}
\item{rasterLegend}{Leave room to the right for 
the legend produced by plotting a Raster object}
\item{\dots}{
Additional arguments passed to \code{legend}, for examaple \code{bty='n'} 
}
}
\details{
\code{scaleBar} produces a scale bar reflecting the distance travelling on a great circle
from the centre of the plot and travelling to the right.  The length of the bar is the width 
of 6 characters times \code{scale.cex}.  

\code{map.new} initiates a plot intended to contain a map covering the extent of \code{x}, 
with no margins.
}
\author{
Patrick Brown
}

\examples{

Npoints = 20
myPoints = SpatialPointsDataFrame(20*cbind(runif(Npoints), runif(Npoints)), 
	data=data.frame(y1=c(NA, rnorm(Npoints-1)), 
	y2=c(sample(1:4, Npoints-1,replace=TRUE), NA)),
	proj4string=CRS("+init=epsg:4326")
	)

\dontrun{
mymap = openmap(myPoints)
}

breaks = c(-100, -1, 1, Inf)
thecol = c('red','orange','blue')


map.new(myPoints)
\dontrun{
plot(mymap,add=TRUE)
}
plot(myPoints,col = as.character(cut(
	myPoints$y1, breaks, thecol
)),add=TRUE)
scaleBar(myPoints, "bottomright",cex=1.25, scale.cex=2)
legendBreaks("topleft", breaks=breaks, col=thecol)

}

