\name{netherlands}
\alias{netherlands}
\alias{nldElev}
\alias{nldTiles}
\alias{nldCities}
\alias{nldCmap}
\alias{meuse}

\docType{data}
\title{
Data from the Netherlands
}
\description{
Elevation data and map tiles for the Netherlands
}
\usage{data("netherlands")}
\format{
\code{nldElev} is a raster of elevation
\code{nltTiles} is a background map
\code{meuse} classic Meuse river data set from the sp package 
\code{nldCities} is a SpatialPointsDataFrame of city locations.
}
\details{
The inclusion of these datasets is intended to allow the package to build 
when an internet connection is not present.
}



\examples{

meuse = terra::unwrap(meuse)
nldTiles = terra::unwrap(nldTiles)
nldCities = terra::unwrap(nldCities)

oldpar=map.new(meuse, buffer=1*1000)
plot(nldTiles,add=TRUE)
points(nldCities, pch=4, col='blue')
text(nldCities,label=nldCities$name, pos=2, col='blue')
points(meuse, pch=15, col=as.integer(meuse$soil))

legend('topleft', fill=1:nlevels(meuse$soil),
	legend=levels(meuse$soil), inset=0.2, bg='white', title='Soil type')
par(oldpar)


}


\keyword{datasets}
