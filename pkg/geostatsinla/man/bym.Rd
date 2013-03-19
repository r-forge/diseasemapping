\name{bym}
\alias{bym}
\alias{bym.SpatialPolygonsDataFrame}
\alias{bym.data.frame}

\title{
Fit the BYM model
}
\description{
Uses inla to fit a Besag, York and Mollie disease mapping model
}
\usage{
bym(data, ...)
	\method{bym}{SpatialPolygonsDataFrame}(data, 
		formula=observed ~ offset(logExpected), 
		priorCI=list(sdSpatial = c(0.01, 10), sdIndep = c(0.01, 10)), 
		family="poisson",
		region.id, adjMat, formula.fitted=formula,...)
	\method{bym}{data.frame}(data, formula, 
		priorCI=list(sdSpatial = c(0.01, 10), sdIndep = c(0.01, 10)), 
		family="poisson",
		region.id, adjMat, formula.fitted=formula,...)
	}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{data}{
The observations and covariates for the model, can be output from 
 \code{\link[diseasemapping]{getSMR}}.
}
\item{formula}{model formula, defaults to intercept-only model suitable for
 output from  \code{\link[diseasemapping]{getSMR}} if \code{data} is a
 SpatialPointsDataFrame. }
\item{priorCI}{prior 95 pct intervals for the standard deviation of the spatial and
independent random effects.
}
\item{family}{distribution of the observations, defaults to \code{poisson}}
\item{region.id}{Variable in \code{data} giving identifiers for the spatial regions.  
If not supplied, row numbers will be used.}
\item{adjMat}{An object of class \code{nb} containing the adjacency matrix.  If not supplied
it will be computed from \code{data}, but is required if \code{data} is a SpatialPolygonDataFrame}
\item{formula.fitted}{formula to use to compute the fitted values, defaults to the model 
formula but may, for example, exclude individual-level covariates.}
  \item{\dots}{
Additional arguments passed to  \code{\link[INLA]{inla}}
}
}
 
\value{
A list containing
\item{inla}{results from the call to \code{\link[INLA]{inla}} }
\item{data}{A data frame or SpatialPolygonsDataFrame containing posterior means
and quantiles of the spatial random effect and fitted values.}
\item{parameters}{summaries of posterior distributions, see \code{\link{glgm}}}
}

\author{
Patrick Brown
}


\seealso{
\code{\link{glgm}},  \code{\link[INLA]{inla}}, \code{\link[diseasemapping]{getSMR}} 
}
\examples{
library(diseasemapping)
data(popdata)
data(casedata)

model = getRates(casedata, popdata, ~age*sex, breaks=seq(0, 90, by=10) )
ontario = getSMR(popdata,model,  casedata,regionCode="CSDUID", regionCodeCases="CSD2006")

res = bym(ontario, region.id="CSDUID", 
formula= observed ~ medIncome + offset(logExpected),
priorCI = list(sdSpatial = c(0.5, 2), sdIndep=c(0.1, 4))
)

res$par$summary

spplot(res$data, "fitted.exp",at=c(0, 0.5, 1, 2,5, 10, 25),
col.regions= terrain.colors(6),xlim=c(-84,-74), ylim=c(41, 48))

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ Bayesian inference }
\keyword{ disease mapping }% __ONLY ONE__ keyword per line
