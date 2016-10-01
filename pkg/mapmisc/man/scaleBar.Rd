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
			cropInset = extent(-180,xmax=180, ymin=-47, ymax=71),
			outer=TRUE, ...) 
}

\arguments{
  \item{crs}{
A \code{CRS} object, proj4 string, or an object from which a projection 
can be extracted with \code{proj4string(crs)}
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
scaleBar(myPoints, "bottomright",cex=1.25, seg.len=2)
temp=legendBreaks("topleft", legend=breaks, col=thecol)

\dontrun{
thedot = insetMap(myPoints, "bottomleft",col='#00000000', lty=0)
points(thedot)
}
}

