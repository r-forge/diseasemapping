\name{glgm}
\alias{glgm}

\title{
Generalized Linear Geostatistical Models
}
\description{
Fits a generalized linear geostatistical model
using   \code{\link[inla]{inla}} }
\usage{
glgm(formula, data, cells, covariates=NULL, priorCI=NULL, ...)
}
 
\arguments{
\item{formula}{Model formula, the spatial random effect should not be supplied but the default 
can be overridden with a 
 {\tt f(space,..)} term. 
}
\item{data}{
An object of class{\tt SpatialPointsDataFrame} containing the data.
}
\item{cells}{Either an integer giving the number of cells in the x direction, or a raster object which 
will be used for the spatial random effect.  If the cells in the raster are not square, the resolution in the y direction 
will be adjusted to make it so.}
\item{covariates}{Either a single raster, a list of rasters or a raster stack containing covariate values used when 
making spatial predictions.  Names of the raster layers or list elements correspond to names in the formula.  If
a covariate is missing from the data object it will be extracted from the rasters.}
}
\item{priorCI}{list with named elements of 0.025 and 0.975 quantiles of prior distributions, or a single 
vector giving the prior quantiles for the range parameter.  List elements can be named: \code{range} for 
the range parameter (not the scale); \code{sd} for the standard deviation (not the variance or precision); 
\code{intercept} for the intercept; and \code{betas} for the regression coefficients.}
}
\details{
This function performs Bayesian inference for generalized linear geostatistical models with INLA.  The Markov random field
approximation on a regular lattice is used for the spatial random effect.
}
\value{
A list with two componnents named \code{raster} and \code{inla}.  \code{inla} contains the results of the call to the
\code{\link[inla]{inla}} function.  \code{raster} is a raster stack with the following layers:
\item{random.mean, random.sd,random.025quant, random.05quant, random.975quant, random.kld}{}
\item{predict.mean, predict.sd,predict.025quant, predict.05quant, predict.975quant, predict.kld}{}
\item{predict.exp}{}
\item{predict.invlogit}{}
}


\seealso{
\code{\link[inla]{inla}}
}
\examples{
data(swissRain)
swissRain$lograin = log(swissRain$rain)
swissFit =  glgm(lograin ~ elevation, swissRain, cells=20,
	covariates=swissAltitude, family="gaussian",
	priorCI=list(sd=c(0.05, 5), range=c(100,1000000),intercept=log(c(1, 1000)), 
	control.mode=list(mode=c(1,1000))

swissFit$inla$summary.fixed 

plot(swissFit$raster[["predict.exp"]]) 



\dontrun{
load(url("http://www.filefactory.com/file/frd1mhownd9/n/CHE_adm0_RData"))
plot(spTransform(gadm, swissRain@proj4string), add=T)
}

}

 

