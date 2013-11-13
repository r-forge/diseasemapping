\name{netherlands}
\alias{netherlands}
\alias{nldElev}
\alias{nldTiles}
\alias{nldCities}
\alias{world}

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
\code{nldCities} is a SpatialPointsDataFrame of city locations.
}
\details{
The inclusion of these datasets is intended to allow the package to build 
when an internet connection is not present.
}

\source{
See examples.  
}

\seealso{\link[sp]{meuse}, \link[raster]{getData}, \link{openmap}}

\examples{

# soil data
library("sp")
data("meuse")
data(meuse)
coordinates(meuse) <- ~x+y
proj4string(meuse) <- CRS("+init=epsg:28992")
meuse$soilFac = factor(meuse$soil, levels=c(1,2,3), 
	labels=c("Calcareous","Non-Calc's","Red Brick"))

\dontrun{
# this is how the data were obtained

# map tiles
nldTiles = openmap(meuse, path="osm-no-labels", zoom=12)

# elevation data
meuseLL = spTransform(meuse, CRS("+init=epsg:4326"))
getData("SRTM", lon=xmin(extent(meuseLL)), 
	lat=ymin(extent(meuseLL)),path=tempdir())
nldElev = raster(paste(tempdir(), "/", "srtm_38_02.tif", sep=""))
nldElev = crop(nldElev, extend(extent(meuseLL), 0.1))
nldElev = projectRaster(nldElev, crs=proj4string(meuse))
nldElev = crop(nldElev, extent(nldTiles))

# cities
nldCities = GNcities(nldTiles, maxRows=25)

# world
world = openmap(extent(-10,30,40,60),path="osm-no-labels")
# save the files where the package builder wants them
save(nldElev, nldTiles, nldCities,world,
file="~/workspace/diseasemapping/pkg/mapmisc/data/netherlands.RData",
compress="xz")
}
data("netherlands")

plot(nldElev)
points(meuse)

map.new(meuse)
plot(nldTiles,add=TRUE)
points(nldCities)
text(nldCities,label=nldCities$name, pos=3)
# insetMap(meuse, "bottomright",map=world)
}


\keyword{datasets}
