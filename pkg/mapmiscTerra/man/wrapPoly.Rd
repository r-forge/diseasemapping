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
\item{crs}{Caracther string representing a  \code{\link[terra]{crs}}.}
}

\value{
A reprojected Spatial object.
}

\seealso{\code{\link[terra]{project}}.}

