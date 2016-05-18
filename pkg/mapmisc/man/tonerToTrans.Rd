\name{tonerToTrans}
\alias{tonerToTrans}
\alias{rgbtToIndex}
\title{
 Convert RGB maps to semi-transparent
}
\description{
Stamen-toner maps are 3-layer RGB rasters, which are converted
to single-layer rasters with indexed colours with whites becoming transparent. 
}
\usage{
tonerToTrans(x,  pattern="(red|green|blue)$", power = 0.5, col='black')
rgbtToIndex(x,  pattern="(red|green|blue|trans)$")
}
\arguments{
  \item{x}{
	A \code{RasterStack} with RGB colours, such as from \code{\link{openmap}} with \code{path='stamen-toner'}}
	\item{pattern}{string passed to \code{\link[base]{grep}} to find RGB layers.}
  \item{power}{
	Values below 1 increase opacity, above 1 increases transparency
}
\item{col}{colour for resulting map}
}

\details{
The difference between these functions is that \code{tonerToTrans} converts white to transparent, whereas \code{rgbtToIndex} uses the transparency layer.  The former is intended for \code{'stamen-toner'} maps.
}

\value{
A \code{RasterLayer} with indexed colours
}

\author{
Patrick Brown
}


\seealso{
\code{\link{openmap}}
}
\examples{

\dontrun{

	rgbMap = openmap(c(0,10), zoom=3, path='stamen-toner')
	names(rgbMap)
	plotRGB(rgbMap)
	
	transMap = tonerToTrans(rgbMap, col='blue')
	names(transMap)
	par(bg='red')
	plot(transMap)

	rgbMap[['stamen.tonerTrans']] = 255-rgbMap[['stamen.tonerRed']]
	rgbtMap = rgbtToIndex(rgbMap)
	plot(rgbtMap)
}
}

