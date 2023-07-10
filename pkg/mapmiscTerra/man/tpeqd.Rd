\name{tpeqd}
\alias{tpeqd}
\alias{tpers}
\title{Two point equidistant and tilted perspective projections}
\description{
Defines map projection
}
\usage{
tpeqd(x, offset=c(0,0), axis='enu')
tpers(x,  hKm = 100*1000, tilt = -10, offset=c(0,0), axis='enu')
}

\arguments{
	\item{x}{A SpatialPoints* object of length 2 or a 
	matrix with two columns. }
	\item{hKm}{Height veiwing the Earth from}
	\item{tilt}{Viewing angle}
\item{offset}{2 coordinates to define the origin}
\item{axis}{defaults to east, north, up.  'swu' would rotateo 90 degrees}
}

\details{
A coordinate reference system is returned
}


\value{
	Caracther string representing a \code{\link[terra:crs]{crs}}.
}	
\references{
\url{https://en.wikipedia.org/wiki/Two-point_equidistant_projection}
\url{https://proj.org/operations/projections/tpers.html}
}
\seealso{\code{\link[terra:crs]{crs}},\code{\link[terra]{project}}}

\examples{

tpeqd(rbind(c(0,0), c(10,50)))

data('netherlands')
nldCities = terra::unwrap(nldCities)

tcrs = tpeqd(x=nldCities[1:2,])

nldT = project(nldCities, tcrs)


map.new(nldT)
text(nldT,labels=nldT$name)
scaleBar(nldT, 'topright')

}

