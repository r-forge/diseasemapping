\name{worldMap}
\alias{worldMap}
\alias{worldmap}
\alias{isohedron}

\docType{data}
\title{
Country boundaries
}
\description{
Country borders from naturalearthdata.com
}
\usage{data("worldMap")}


\source{
\url{https://www.naturalearthdata.com/downloads/110m-cultural-vectors/}}


\examples{

# soil data

data("worldMap")
worldMap = terra::unwrap(worldMap)
oldpar=map.new(worldMap)
plot(worldMap, border='red', lwd=3, add=TRUE)
plot(worldMap[worldMap$NAME ==  'Brazil',],
	add=TRUE, col='green')
par(oldpar)

}


\keyword{datasets}
