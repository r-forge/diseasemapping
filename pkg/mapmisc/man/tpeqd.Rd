\name{tpeqd}
\alias{tpeqd}
\title{Two point equidistant projections}
\description{
Defines map projection
}
\usage{
tpeqd(x, offset=c(0,0))
}

\arguments{
	\item{x}{A SpatialPoints* object of length 2 or a 
	matrix with two columns. }
\item{offset}{2 coordinates to define the origin}
}

\details{
A coordinate reference system is returned
}


\value{
	An object of class \code{\link[sp]{CRS}}.
}	

\references{
\url{http://en.wikipedia.org/wiki/Two-point_equidistant_projection}
}
\seealso{\code{\link[sp]{CRS}},\code{\link[sp]{spTransform}}}

\examples{

tpeqd(rbind(c(0,0), c(10,50)))

data('netherlands')

tcrs = tpeqd(nldCities[1:2,])
tcrs


nldT = spTransform(nldCities, tcrs)
projection(nldT)

map.new(nldT)
text(nldT,labels=nldT$name)
scaleBar(nldT, 'topright')



}