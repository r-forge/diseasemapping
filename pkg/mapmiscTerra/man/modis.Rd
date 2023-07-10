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
getModisTiles(x, tiles)
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

myPointLL = vect(cbind(c(5:6),10:11), crs = crsLL)


  getModisTiles(myPointLL)

  getModisTiles(myPointLL, degreeRaster)



modisUrl = 	'https://e4ftl01.cr.usgs.gov/MOTA/MCD12Q1.006/2002.01.01/'

\dontrun{
theFiles = 	grep(
		paste0("(",paste(getModisTiles(myPointLL)[,'tile'], collapse='|'), ").*.hdf$"),
	XML::getHTMLLinks(RCurl::getURL(
	    modisUrl,ftp.use.epsv=FALSE,
	    dirlistonly = TRUE)), value=TRUE)
paste0(modisUrl, theFiles)
}
}
