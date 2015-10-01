\name{wrapPoly}
\alias{wrapPoly}

\title{Reproject polygons with wrapping}
\description{
Reprojects a SpatialPolygons object to a projection with longitude wrapping other than 180 degreess
}
\usage{
wrapPoly(x,crs)
	}

\arguments{
\item{x}{A Spatial object}
\item{crs}{	An object of class \code{\link[sp]{CRS}}.}
}

\value{
A reprojected Spatial object.
}

\seealso{\code{\link[sp]{spTransform}}.}

