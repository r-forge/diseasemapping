\name{geostatData}
\alias{geostatData}
\alias{geostatData.default}
\alias{geostatData.SpatRaster}
\alias{geostatData.SpatVector}

\title{Prepare observations and covariates for a geostatistical model}
\description{
Extracts covariates at observation locations, aligns covariate rasters to a
prediction grid, and returns the data structures used by \code{\link{glgm}}
and \code{\link{lgm}}.
}
\usage{
geostatData(formula, data, grid, covariates, buffer = 0)
\method{geostatData}{SpatVector}(formula, data, grid, covariates, buffer = 0)
\method{geostatData}{SpatRaster}(formula, data, grid = data, covariates = NULL, buffer = 0)
\method{geostatData}{default}(formula, data, ...)
}

\arguments{
\item{formula}{Model formula.  Names on the right-hand side are taken from
\code{data} when present and otherwise extracted from \code{covariates}.}
\item{data}{A \code{SpatVector} of point observations, or a \code{SpatRaster}
of observations on a lattice.}
\item{grid}{A \code{SpatRaster} giving the prediction grid (and the lattice
for the spatial random effect).  For the \code{SpatRaster} method this
defaults to \code{data}.}
\item{covariates}{A \code{SpatRaster}, a multi-layer \code{SpatRaster}, or a
list of \code{SpatRaster}s.  Layer or list names must match names in
\code{formula}.}
\item{buffer}{Extra space padded around \code{grid}, passed to the same
argument of \code{\link{squareRaster}}.}
\item{\dots}{Ignored.  The default method reports an error if \code{data} is
not a \code{SpatVector} or \code{SpatRaster}.}
}

\value{
A list with
\item{data}{For \code{SpatVector} input, a \code{SpatVector} with extracted
covariates and a \code{space} cell identifier.  For \code{SpatRaster} input,
a \code{data.frame} of cell values.}
\item{grid}{The (possibly buffered) square \code{SpatRaster} used as the
model lattice.}
\item{covariates}{A \code{data.frame} of covariate values on \code{grid}.}
\item{formula}{Returned only by the \code{SpatRaster} method.  Equal to the
supplied formula, or a rewritten formula when log-offsets are aggregated.}
}

\seealso{
\code{\link{glgm}}, \code{\link{lgm}}, \code{\link{squareRaster}}
}

\examples{
data("loaloa")
loaloa = unwrap(loaloa)
elevationLoa = unwrap(elevationLoa)
eviLoa = unwrap(eviLoa)

res = geostatData(
    y ~ elev + evi,
    data = loaloa,
    covariates = list(elev = elevationLoa, evi = eviLoa),
    grid = squareRaster(loaloa, cells = 20, buffer = 1e4)
)
names(res)
names(res$data)
res$grid
}
