\name{modis}
\alias{modis}
\alias{crsModis}
\alias{modisRaster}
\alias{getModisTiles}
\alias{degreeRaster}
\title{MODIS tiles and projection}
\description{
Raster containing MODIS tile ID's
}
\usage{
getModisTiles(x, tiles = mapmisc::modisRaster)
crsModis
modisRaster
degreeRaster
}
\arguments{
\item{x}{A spatial object which modis tiles will cover.}
\item{tiles}{A raster with modis (or other) tiles.}
}

\details{
	\code{getModisTiles} returns a matrix with modis tiles.

	\code{modisRaster} shows horizontal and vertical tile names for downloading data from MODIS
	at, for example, \url{http://e4ftl01.cr.usgs.gov/MOTA}.
	
	\code{degreeRaster} shows horizontal and vertical tiles in long-lat, for downloading elevation.
}



\references{
\url{https://nsidc.org/data/docs/daac/mod10_modis_snow/landgrid.html},
\url{https://spatialreference.org/ref/sr-org/modis-sinusoidal/}
}


\examples{

crsModis

myPointLL = SpatialPoints(cbind(c(5:6),10:11), proj4string = crsLL)

if(requireNamespace("rgdal", quietly=TRUE)) {
  getModisTiles(myPointLL)

  getModisTiles(myPointLL, mapmisc::degreeRaster)
}

modisUrl = 	'ftp://ladsweb.nascom.nasa.gov/allData/5/MCD12Q1/2002/001/'

\dontrun{
paste(modisUrl,
  grep(
    paste(getModisTiles(myPointLL)[,'tile'], collapse='|'),
    unlist(strsplit(RCurl::getURL(
	    modisUrl,ftp.use.epsv=TRUE,
	    dirlistonly = TRUE), '\\n')), 
	  value=TRUE), sep='')
}
}
