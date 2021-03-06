\name{GNcities}
\alias{GNcities}
\alias{GNsearch}
\title{Retrieve city names and locations}
\description{
 This function uses the geonames package to  provide city names and locations
 from  www.geonames.org.
  }
\usage{
	GNcities(north, east, south, west, lang = "en", maxRows = 10, buffer=0)
	GNsearch(..., crs=crsLL)
	}

\arguments{
	\item{north}{
	A bounding box or  
	SpatialPoints or SpatialPolygons or 
	Extent or  Raster object, or a decimal degree of longitude.}
	\item{east,south,west}{If \code{north} is numeric, decimal degree bounding box.}
	\item{lang}{Language for internationalised returned text}
	\item{maxRows}{Limit on returned rows}
	\item{buffer}{passed to code{\link[raster]{extend}}}
	\item{...}{Various search arguments}
	\item{crs}{projection for the output}
}

\seealso{\code{\link[geonames]{GNcities}}, \code{\link[geonames]{GNsearch}}}


\value{A SpatialPointsDataFrame with the sampe projection  \code{north} if it exists, otherwise in 
long-lat.}


\examples{

myraster = raster(matrix(0,10,10),xmn=8,xmx=18,ymn=0,ymx=10, 
	crs=mapmisc::crsLL)
values(myraster) = seq(0,1,len=ncell(myraster))
myPoints = SpatialPoints(myraster, proj4string=CRS(proj4string(myraster)))[
	seq(1,ncell(myraster),len=5)]

\dontrun{

if (requireNamespace("geonames", quietly = TRUE)) { 
cities=GNcities(myPoints, max=5)
mytiles = openmap(myraster)

map.new(cities)
plot(mytiles, add=TRUE)
points(cities, col='red')
text(cities, labels=cities$name, col='red',pos=4)


cities=GNcities(myraster, max=5)

map.new(cities)
plot(mytiles, add=TRUE)
points(cities, col='red')
text(cities, labels=cities$name, col='red',pos=4)

mapmisc::GNsearch(q="Toronto Ontario")
}
}


}