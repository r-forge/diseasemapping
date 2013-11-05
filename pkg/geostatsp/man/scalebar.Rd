\name{scalebar}
\alias{scalebar}
\alias{map.new}
\alias{legendBreaks}

\title{
Map, scale bar, and legend
}
\description{
Utilities for plotting a map, adding a scale bar and north arrow, 
and adding a legend of colour scales.
}
\usage{
scalebar(crs, pos = "bottomright", scale.cex = 1, ...)
legendBreaks(pos, col, breaks, ...)
map.new(x)
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
\item{col}{A vector of colours}
\item{breaks}{A vector of breaks, should be have length \code{1+length(col)}}
\item{x}{
A spatial object from which an extent can be extracted.
}
\item{\dots}{
Additional arguments passed to \code{legend}, for examaple \code{bty='n'} 
}
}
\details{
\code{scalebar} produces a scale bar reflecting the distance travelling on a great circle
from the centre of the plot and travelling to the right.  The length of the bar is the width 
of 6 characters times \code{scale.cex}.  

\code{legendBreaks} produces a legend for 'z-axis' colour scales.

\code{map.new} initiates a plot intended to contain a map covering the extent of \code{x}, 
with no margins.
}
\author{
Patrick Brown
}

\examples{

data("swissRain")
map.new(swissRain)

mybreaks = c(0,  1000, 2000, 3000, 4500)
mycol = c("lightgrey","yellow","orange","red")

plot(swissAltitude, col=mycol,breaks=mybreaks,add=TRUE,legend=FALSE)
points(swissRain)
plot(swissBorder, add=TRUE)
scalebar(swissRain, "bottomright",cex=1.25, scale.cex=2)
legendBreaks("topleft", col=mycol, breaks=mybreaks,title='elevation')


}

