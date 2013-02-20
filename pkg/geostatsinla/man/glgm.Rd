\name{glgm}
\alias{glgm}
\alias{lgcp}
\alias{geostatsinla-package}

\title{
Generalized Linear Geostatistical Models
}
\description{
Fits a generalized linear geostatistical model or a log-Gaussian Cox process
using   \code{\link[inla]{inla}} }
\usage{
glgm(data,  cells, covariates=NULL, formula=NULL, priorCI=NULL, maternRoughness=1, mesh=F,...)
lgcp(data,  cells, covariates=NULL, formula=NULL, priorCI=NULL, maternRoughness=1, mesh=F,...)
}
 
\arguments{
\item{data}{
An object of class {\tt SpatialPointsDataFrame} containing the data.
}
\item{cells}{Either an integer giving the number of cells in the x direction, or a raster object which 
will be used for the spatial random effect.  If the cells in the raster are not square, the resolution in the y direction 
will be adjusted to make it so.}
\item{covariates}{Either a single raster, a list of rasters or a raster stack containing covariate values used when 
making spatial predictions.  Names of the raster layers or list elements correspond to names in the formula.  If
a covariate is missing from the data object it will be extracted from the rasters.  Defaults to \code{NULL} for an 
intercept-only model.
}
\item{formula}{Model formula, defaults to a linear combination of each of the layers in the \code{covariates} object.   
The spatial random effect should not be supplied but the default 
can be overridden with a 
 {\tt f(space,..)} term. For \code{glgm} the response variable defaults to the first variable in the \code{data} object, and
 \code{formula} can be an integer or character string specifying the response variable.  For \code{lgcp}, the formula
 should be one-sided.
}
\item{priorCI}{list with named elements of 0.025 and 0.975 quantiles of prior distributions, or a single 
vector giving the prior quantiles for the range parameter.  List elements can be named: \code{range} for 
the range parameter (not the scale); \code{sd} for the standard deviation (not the variance or precision); 
\code{intercept} for the intercept; and \code{betas} for the regression coefficients.}
\item{maternRougness}{roughness parameter for the Matern correlation function, must be 1 or 2.}
\item{mesh}{Currently unimplemented options for using a mesh rather than a grid for the Markov random field approximation.}
\item{...}{Additional options passed to  \code{\link[inla]{inla}}}
}
\details{
This function performs Bayesian inference for generalized linear geostatistical models with INLA.  The Markov random field
approximation on a regular lattice is used for the spatial random effect.
}
\value{
A list with two components named  \code{inla}, \code{raster}, and \code{parameters}.  \code{inla} contains the results of the call to the
\code{\link[inla]{inla}} function.  \code{raster} is a raster stack with the following layers:
\item{random.mean, random.sd,random.025quant, random.05quant, random.975quant, random.kld}{}
\item{predict.mean, predict.sd,predict.025quant, predict.05quant, predict.975quant, predict.kld}{}
\item{predict.exp}{}
\item{predict.invlogit}{Only supplied if a binomial response variable was used.}
\code{parameters} contains a list with element
\item{summary}{}
and \code{range} and \code{sd} elements containing, 
for the range and standard deviation parameters respectively,
\item{prior}{}
\item{posterior}{}
}

\seealso{
\code{\link[inla]{inla}}
}
\examples{
# geostatistical model for the swiss rainfall data
library(geostatsp)
data(swissRain)
swissRain$lograin = log(swissRain$rain)
library(INLA)
swissFit =  glgm(swissRain, cells=20, formula="lograin",
	covariates=swissAltitude, family="gaussian",
	priorCI=list(sd=c(0.5, 1.5), range=c(150000,250000),intercept=log(c(1, 1000))), 
	control.mode=list(theta=c(0.5,0.1,1.5)),
	control.family=list(hyper=list(prec=list(prior="loggamma", param=c(20, 20))))
	)

swissFit$inla$summary.fixed

swissFit$parameters$summary

plot(swissFit$raster[["predict.exp"]]) 

\dontrun{
load(url("http://www.filefactory.com/file/frd1mhownd9/n/CHE_adm0_RData"))
plot(spTransform(gadm, swissRain@proj4string), add=T)
}

# a log-Gaussian Cox process example

}

 
data=swissRain; cells=20; formula="lograin"
	covariates=swissAltitude; family="gaussian"
	priorCI=list(sd=c(0.9, 1),  range=c(190000,200000),
	intercept=log(c(1, 1000))) 
	control.mode=list(mode=c(1,1000))
		control.family=list(hyper=list(prec=list(prior="loggamma", param=c(20, 20))))
	maternRoughness=2