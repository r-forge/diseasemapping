\name{tonerToTrans}
\alias{tonerToTrans}
\title{
 Convert RGB maps to semi-transparent
}
\description{
Stamen-toner maps are 3-layer RGB rasters, which are converted
to single-layer rasters with indexed colours with whites becoming transparent. 
}
\usage{
tonerToTrans(x,  pattern="(red|green|blue)$", power = 0.5, col='black',  threshold=Inf)
}
\arguments{
  \item{x}{
	A \code{RasterStack} with RGB colours, such as from \code{\link{openmap}} with \code{path='stamen-toner'}}
	\item{pattern}{string passed to \code{\link[base]{grep}} to find RGB layers.}
  \item{power}{
	Values below 1 increase opacity, above 1 increases transparency
}
\item{col}{colour for resulting map}
\item{threshold}{colours above this value are transparent}
}

\details{
Intended for \code{'stamen-toner'} maps.
}

\value{
A \code{SpatRast} with indexed colours
}

\author{
Patrick Brown
}


\seealso{
\code{\link{openmap}}
}
\examples{

	data('worldMap');worldMap = unwrap(worldMap)

	origMap = openmap(worldMap[worldMap$NAME == 'Sierra Leone',],
	 path='stamen-toner',  
		buffer=100*1000
	)
		# path = 'https://tiles.stadiamaps.com/tiles/stamen_toner/')
	#	suffix = '@2x.png')
	oldpar=	map.new(origMap, bg='green')
	plot(origMap, add=TRUE)
	
	transMap = tonerToTrans(origMap)
	names(transMap)
  map.new(transMap, bg='green')
	plot(transMap, add=TRUE)
	par(oldpar)


}

