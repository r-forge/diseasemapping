\name{modis}
\alias{modis}
\alias{crsModis}
\alias{modisRaster}
\title{MODIS tiles and projection}
\description{
Raster containing MODIS tile ID's
}
\usage{
crsModis
modisRaster
}


\details{
	\code{modisRaster} shows horizontal and vertical tile names for downloading data from MODIS
	at, for example, \url{http://e4ftl01.cr.usgs.gov/MOTA}.
}



\references{
\url{https://modis-land.gsfc.nasa.gov/MODLAND_grid.html},
\url{http://spatialreference.org/ref/sr-org/6842}
}


\examples{

crsModis

	myPointLL = SpatialPoints(
		cbind(5,10), proj4string = crsLL
	)

if(require('rgdal', quietly=TRUE)) {
	myPointModis = spTransform(myPointLL, crsModis)
} else {
	myPointModis = SpatialPoints(
		cbind(547500, 1111200), 
		proj4string = crsModis
	)
}

factorValues(modisRaster,
	extract(modisRaster, myPointModis)
)

paste(
  'ftp://ladsweb.nascom.nasa.gov/allData/5/MCD12Q1/2002/001/MCD12Q1.A2002001.',
  factorValues(modisRaster,
    values(raster::crop(
	    modisRaster, extend(extent(myPointModis), 1500*1000)
    )))[,'tile'],
  '.005.*', sep='')    

}