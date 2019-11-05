\name{worldMap}
\alias{worldMap}
\alias{worldmap}

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
library("mapmisc")

data("worldMap")
myMap = openmap(worldMap, zoom=1, fact=2)
map.new(worldMap)
plot(myMap, add=TRUE)
plot(worldMap, border='red', lwd=3, add=TRUE)
plot(worldMap[worldMap$NAME ==  'Brazil',],
	add=TRUE, col='green')

}


\keyword{datasets}
