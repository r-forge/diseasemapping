\name{crsMerc}
\alias{crsMerc}
\alias{crsCanada}
\alias{crsLL}
\alias{bboxLL}
\alias{bboxLLsafe}
\alias{extentMerc}


\title{Some coordinate reference systems and bounding boxes}
\description{
Defines CRS's for the several map projections.
}
\usage{
crsMerc 
crsLL 
crsCanada 
extentMerc
bboxLLsafe
bboxLL
}


\format{
\code{crsMerc} spherical Mercator projection used by web mapping services, epsg:3857
\code{crsLL} long-lat, epsg:4326 
\code{crsCanada} customized oblique mercator for Canada
\code{bboxLL} polygon of bounding box of long-lat, -180 to 180, -90 to 90
\code{bboxLLsafe} as bboxLL, but slightly away from the edges
\code{extentMerc} extent of spherical mercator projections
}

\details{
	these objects are used internally and may be of interest to the user
}

\value{objects of class \code{\link[terra:crs]{crs}} or numeric vectors.}

\references{
\url{https://en.wikipedia.org/wiki/Web_Mercator},
\url{https://spatialreference.org/ref/epsg/4326/}
}
\seealso{\code{\link[terra:crs]{crs}}}

\examples{
terra::crs(crsMerc, proj=TRUE)
terra::crs(crsLL, proj=TRUE)
terra::crs(crsCanada, proj=TRUE)
terra::ext(extentMerc)

bboxLLsafe = terra::unwrap(bboxLLsafe)
plot(bboxLLsafe)
plot(terra::project(bboxLLsafe, crsMerc))

}