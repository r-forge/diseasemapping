\name{tonerToTrans}
\alias{tonerToTrans}
\title{
 Convert stamen-toner maps to transparent
}
\description{
Stamen-toner maps are 3-layer RGB rasters, which are converted
to single-layer rasters with indexed colours with whites becoming transparent. 
}
\usage{
tonerToTrans(x, power = 0.5)
}
\arguments{
  \item{x}{
	A \code{RasterStack} from \code{\link{openmap}} with \code{path='stamen-toner'}
}
  \item{power}{
	Values below 1 increase opacity, above 1 increases transparency
}
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
	
	transMap = tonerToTrans(rgbMap)
	names(transMap)
	par(bg='red')
	plot(transMap)

}

