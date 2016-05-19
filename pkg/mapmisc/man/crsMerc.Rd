\name{crsMerc}
\alias{crsMerc}
\alias{crsLL}
\title{Spherical Mercator projection}
\description{
Defines CRS's for the Spherical Mercator and long-lat projections.
}
\usage{
crsMerc
crsLL
}


\details{
CRS objects for epsg:4326 (long-lat) and the spherical Mercator projection used by web mapping services.
Using epsg codes requires the \code{rgdal} package to be installed, and 
\code{crsLL} is intended as a replacement for \code{CRS("+init=epsg:4326")} when 
\code{rgdal} is not guaranteed to be available.  
}


\value{
	Objects of class \code{\link[sp]{CRS}}.
}	

\references{
\url{https://en.wikipedia.org/wiki/Web_Mercator},
\url{http://spatialreference.org/ref/epsg/4326/}
}
\seealso{\code{\link[sp]{CRS}}}

\examples{

crsMerc
try(CRS("+init=epsg:3857"), silent=TRUE)

crsLL
try(CRS("+init=epsg:4326"), silent=TRUE)

if(require('rgdal', quietly=TRUE)){
	xOrig = SpatialPoints(cbind(-40,30), proj4string=crsLL)
	xMerc = spTransform(xOrig, crsMerc)
	xLL = spTransform(xMerc, crsLL)
	
	coordinates(xOrig)
	coordinates(xMerc)
	coordinates(xLL)
}
	
}