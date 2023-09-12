\name{modis}
\alias{modis}
\alias{crsModis}
\alias{getModisRaster}
\alias{getModisTiles}
\alias{getDegreeRaster}
\title{MODIS tiles and projection}
\description{
Raster containing MODIS tile ID's
}
\usage{
getModisTiles(x, tiles)
crsModis
getModisRaster()
getDegreeRaster()
}
\arguments{
\item{x}{A spatial object which modis tiles will cover.}
\item{tiles}{A raster with modis (or other) tiles.}
}

\details{
	Provides information on tiles which can be downloaded from MODIS.
}



\references{
\url{https://modis-land.gsfc.nasa.gov/MODLAND_grid.html},
\url{https://spatialreference.org/ref/sr-org/modis-sinusoidal/}
}


\value{
	\code{getModisTiles} returns a matrix with modis tiles.

	\code{getModisRaster} shows horizontal and vertical tile names for downloading data from MODIS.
	
	\code{getDegreeRaster} shows horizontal and vertical tiles in long-lat, for downloading elevation.
}

\examples{

crsModis

myPointLL = vect(cbind(c(5:6),10:11), crs = crsLL)


getModisTiles(myPointLL)

getModisTiles(myPointLL, getDegreeRaster())



modisUrl = 	'https://e4ftl01.cr.usgs.gov/MOTA/MCD12Q1.061/2002.01.01/'
desiredTiles = paste0("(",
	paste(getModisTiles(myPointLL, getModisRaster())[,'tile'], collapse='|'), 
	").*.hdf$")

if(requireNamespace("RCurl", quietly=TRUE) & requireNamespace("XML", quietly=TRUE)) {
theFiles = 	grep(desiredTiles,
	XML::getHTMLLinks(RCurl::getURL(
	    modisUrl,ftp.use.epsv=FALSE,
	    dirlistonly = TRUE)), value=TRUE)
paste0(modisUrl, theFiles)
}

}
