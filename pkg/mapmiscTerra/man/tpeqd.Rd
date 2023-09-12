\name{tpeqd}
\alias{tpeqd}
\alias{tpers}
\title{Two point equidistant and tilted perspective projections}
\description{
Defines map projection
}
\usage{
tpeqd(x, offset=c(0,0), axis='enu')
tpers(x,  hKm = 100*1000, tilt = -10,   azi, offset=c(0,0), axis='enu')
}

\arguments{
	\item{x}{A SpatialPoints* object of length 2 or a 
	matrix with two columns. }
	\item{hKm}{Height veiwing the Earth from}
	\item{tilt}{Viewing angle}
	\item{azi}{Azimuth, defaults to direction of first two points in x}
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

data('worldMap');worldMap=unwrap(worldMap)

thepoints = vect(rbind(cbind(150, -40), cbind(-70,-40)), crs=crsLL)
crsOne = tpeqd(thepoints)
worldMapTrans = wrapPoly(worldMap, crsOne)

oldpar=map.new(crsOne, col='lightblue')
plot(worldMapTrans, add=TRUE, col='grey')
points(project(thepoints, crsOne), col='red')
gridlinesWrap(crsOne, col='orange')



thepoints = vect(rbind(cbind(-40, 65), cbind(139,35)), crs=crsLL)
crsTwo = tpeqd(thepoints)


map.new(crsTwo, col='lightblue')
plot(wrapPoly(worldMap, crsTwo), add=TRUE, col='grey')
points(project(thepoints, crsTwo), col='red')
gridlinesWrap(crsTwo, col='orange')

par(oldpar)
}


