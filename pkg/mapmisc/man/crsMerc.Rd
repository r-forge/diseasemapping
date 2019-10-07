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
	Objects of class \code{\link[sp:CRS-class]{CRS}}.
}	

\references{
\url{https://en.wikipedia.org/wiki/Web_Mercator},
\url{https://spatialreference.org/ref/epsg/4326/}
}
\seealso{\code{\link[sp:CRS-class]{CRS}}, \code{\link[rgdal:showWKT]{showEPSG}}}

\examples{

crsMerc
if(require('rgdal', quietly=TRUE))
	CRS("+init=epsg:3857")

crsLL
if(require('rgdal', quietly=TRUE)) {
	CRS("+init=epsg:4326")

	rgdal::showEPSG("+proj=longlat")
	rgdal::showEPSG(as.character(crsLL))
}


}