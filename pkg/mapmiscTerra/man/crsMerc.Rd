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
}


\value{
	Character string suitable for \code{\link[terra:crs]{crs}}.
}	

\references{
\url{https://en.wikipedia.org/wiki/Web_Mercator},
\url{https://spatialreference.org/ref/epsg/4326/}
}
\seealso{\code{\link[terra:crs]{crs}}}

\examples{

crsMerc

crsLL

}