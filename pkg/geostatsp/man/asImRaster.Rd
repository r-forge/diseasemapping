\name{asImRaster}
\alias{asImRaster}


\title{Convert a raster to an im object}
\description{
Conversion between rasters and spatstat's im objects
}
\usage{
	asImRaster(X, ...)
}

\arguments{
	\item{X}{
A RasterLayer}
\item{...}{ additional arguments}
}

\details{
This function is not registered as an S3 method as doing so would require the spatstat package to be a dependency of geostatsp.
}

\examples{


if(require('spatstat',quietly=TRUE)){	
	# numeric raster
	myraster = raster(matrix(100+1:5,10,10),
		xmn=1,xmx=6,ymn=-1,ymx=4)

	myIm = asImRaster(myraster)
	
	# factor raster
	levels(myraster) = data.frame(
		ID = sort(unique(values(myraster))), 
		labels=letters[seq(1, length(unique(values(myraster))))]
	)

	myImF = asImRaster(myraster)

	myImF$v
	levels(myraster)
	
	# loaloa land types
	data('loaloa')
	ltLoaIm = asImRaster(ltLoa)
	ltLoaIm$v[1:4,1:5]
}
}


