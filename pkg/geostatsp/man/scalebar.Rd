\name{scalebar}
\alias{scalebar}

\title{
Add scale bar
}
\description{
Add a scale bar and north arrow to the current plot
}
\usage{
scalebar(crs, pos = "bottomright", scale.cex = 1, ...)
}

\arguments{
  \item{crs}{
A \code{CRS} object, proj4 string, or an object from which a projection 
can be extracted with \code{proj4string(crs)}
}
  \item{pos}{
Position of the scale bar, as specified in the \code{legend} function.
}
  \item{scale.cex}{
Scaling factor for the scale bar.
}
  \item{\dots}{
Additional arguments passed to \code{legend}, for examaple \code{bty='n'} 
}
}
\details{
The scale bar is the distance travelling on a great circle
from the centre of the plot and the mid-point of the right border, 
and roughly \code{scale.cex *0.01} times the distance between these two points.
}
\author{
Patrick Brown
}

\examples{

data("swissRain")
plot(swissAltitude, main="elevation", bigplot=c(0,1,0,1))
points(swissRain)
plot(swissBorder, add=TRUE)
scalebar(swissRain)
scalebar(proj4string(swissBorder), "topleft")


}

