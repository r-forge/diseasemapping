\name{wrapPoly}
\alias{wrapPoly}
\alias{llCropBox}

\title{Reproject polygons with wrapping}
\description{
Reprojects a SpatialPolygons object to a projection with longitude wrapping other than 180 degreess
}
\usage{
wrapPoly(x,crs, buffer.width = 100*1000)
llCropBox(crs, 
  buffer.width=50*1000, densify.interval = 25*1000, 
  crop.distance = 2.1e7, crop.poles = FALSE, crop.leftright=FALSE,
  remove.holes=TRUE, cycles = 2, ellipse=NULL)
}

\arguments{
\item{x}{A Spatial object}
\item{crs}{Caracther string representing a  \code{\link[terra]{crs}}.}
\item{buffer.width}{buffer to add to points on border when cropping poloygons, defaults to 100km}
\item{densify.interval}{interval when densifying}
\item{ crop.distance}{crop coordinates larger than this value}
\item{crop.poles}{remove areas near the poles}
\item{crop.leftright}{remove points near 180 longitute line}
\item{remove.holes}{fill holes in the crop region}
\item{cycles}{iterations adding denser points}
\item{ellipse}{boundary of the world in crs coordinates}
}

\value{
A reprojected Spatial object.
}

\seealso{\code{\link[terra]{project}}, examples in \code{\link{tpeqd}}}

