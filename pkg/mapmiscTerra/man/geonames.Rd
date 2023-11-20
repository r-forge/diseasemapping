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
	\item{buffer}{passed to \code{\link[terra]{extend}}}
	\item{...}{Various search arguments}
	\item{crs}{projection for the output}
}

\seealso{\code{\link[geonames]{GNcities}}, \code{\link[geonames]{GNsearch}}}


\value{A SpatialPointsDataFrame with the sampe projection  \code{north} if it exists, otherwise in 
long-lat.}


\examples{

\dontrun{
	GNsearch(q="Toronto Ontario", maxRows = 3)
}


library('terra')
myraster = rast(
	matrix(1:100,10,10),
	extent=ext(8,18,0,10), crs=crsLL)

options(geonamesUsername="myusernamehere") 
if(file.exists("~/geonamesUsername.R")) source("~/geonamesUsername.R")

if(requireNamespace("geonames", quietly = TRUE)) { 

cities=try(GNcities(myraster, max=5), silent=TRUE)
mytiles = openmap(myraster, zoom=5, buffer=1)

oldpar=map.new(mytiles)
plot(mytiles, add=TRUE)
if(!all(class(cities)=='try-error')) {
points(cities, col='red')
text(cities, labels=cities$name, col='red',pos=4)
}


par(oldpar)

}



}