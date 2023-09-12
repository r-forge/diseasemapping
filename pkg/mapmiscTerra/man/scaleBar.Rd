\name{scaleBar}
\alias{scaleBar}
\alias{insetMap}

\title{
Scale bar and inset map
}
\description{
Utilities for plotting a map, adding a scale bar and north arrow, 
and adding a legend of colour scales.
}
\usage{
scaleBar(crs, pos = "bottomright",   
	cex=1,
    pt.cex = 1.1*cex,
    seg.len=5*cex,
    title.cex=cex,
    outer=TRUE,...)
insetMap(crs, pos="bottomright",map="osm",zoom=0, 
		width=max(c(0.2, 1-par('plt')[2])), 
		col="#FF000090", borderMap=NULL,
			cropInset = terra::ext(-180, 180,  -47, 71),
			outer=TRUE, inset = c(0.1, 0.1), ...) 
}

\arguments{
  \item{crs}{
A character string from which a projection 
can be extracted with \code{terra::crs}
}
  \item{pos}{
Position, as specified in the \code{legend} function.
}
\item{cex}{scaling factor for the legend}
\item{pt.cex}{
Scaling factor north arrow (can be zero).
}
\item{seg.len}{approximate length (in character units) of the scale bar.  can be zero.}
\item{title.cex}{scaling for the distance text}
\item{outer}{If \code{TRUE}, put bar or map in the margin of the plot}
\item{map}{Either a Raster for the inset map or a string 
passed to \code{\link{openmap}}'s path argument}
\item{zoom}{Zoom level if retrieving inset map from \code{\link{openmap}}}
\item{width}{Width of the inset map, as a fraction of the plot window}
\item{col}{Colour for shaded region of inset map}
\item{borderMap}{border style for the inset map (passed to \code{\link[graphics]{polygon}})}
\item{cropInset}{Crop the insert map to this extent}
\item{inset}{how far from the border to put the inset map}
\item{\dots}{
Additional arguments passed to \code{\link[graphics]{legend}} for \code{scaleBar} 
or \code{\link[graphics]{polygon}} (for \code{insetMap}).
}
}
\details{
\code{scaleBar} produces a scale bar reflecting the distance travelling on a great circle
from the centre of the plot and travelling to the right.  The length of the bar is the width 
of 6 characters times \code{scale.cex}.  

}

\value{A list containig coordinates of the elements of the scale bar.}

\author{
Patrick Brown
}


\examples{

Npoints = 20
set.seed(0)
myPoints = vect(
	cbind(runif(Npoints)-0.1, 51+runif(Npoints)), 
	atts=data.frame(
		y1=c(NA, rnorm(Npoints-1)), 
	    y2=c(sample(0:5, Npoints-1,replace=TRUE), NA)
	),
	crs=crsLL)



breaks = c(-100, -1, 1, Inf)
thecol = c('red','orange','blue')


oldpar = map.new(myPoints)
plot(myPoints,col = as.character(cut(
	myPoints$y1, breaks, thecol
)),add=TRUE)
scaleBar(myPoints, "bottomright",cex=1.25, seg.len=2)
legendBreaks("topleft", legend=breaks, col=thecol)

thedot = insetMap(crs=myPoints, pos="bottomleft",
	 col='#00000000', lty=0, outer=FALSE, width=0.5)
points(thedot)

par(oldpar)
}

