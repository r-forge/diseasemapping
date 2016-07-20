\name{lgcpToGeostatsp}
\alias{lgcpToGeostatsp}
\title{Convert lgcp results to geostatsp format}
\usage{
lgcpToGeostatsp(x, dirname = x$gridfunction$dirname)
}
\arguments{
\item{x}{result from \code{\link[lgcp]{lgcp}}}

\item{dirname}{folder with detailed results}
}
\description{
Results from \code{\link[lgcp]{lgcp}} are converted 
to a format comparable to that produced by \code{\link[geostatsp]{lgcp}}
from the \code{geostatsp} package
}

