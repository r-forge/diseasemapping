\name{lglm}
\alias{lglm}
 

\title{
 Linear Geostatistical Models
}
\description{
Fits a   linear geostatistical model using Maximum Likelihood Estimation }
\usage{
lgm(data,  cells, covariates=NULL, formula=NULL, mesh=F,...)
}
 
\arguments{
\item{data}{
An object of class \code{ SpatialPointsDataFrame} containing the data.
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
The spatial random effect should not be supplied. The response variable defaults to the first variable in the \code{data} object, and
 \code{formula} can be an integer or character string specifying the response variable.  For \code{lgcp}, the formula
 should be one-sided.
}
\item{...}{Additional options passed to  \code{\link{likfit}}}
}
\details{

}
\value{
A list with two components named  \code{raster}, and \code{parameters}.  \code{inla} contains the results of the call to the
 \code{raster} is a raster stack with the following layers:
\item{random.mean, random.sd,random.X0.025quant, random.X0.5quant, random.X0.975quant, random.kld}{}
\item{predict.mean, predict.sd,predict.X0.025quant, predict.X0.5quant, predict.X0.975quant, predict.kld}{}
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
\code{\link{likfit}}, \code{\link{krige}}

}
\examples{

data(swissRain)
library(INLA)
swissFit =  lgm(swissRain, cells=40, formula="lograin",
	covariates=swissAltitude, lambda=0.5, fix.lambda=T, 
	fix.psiR=F, fix.psiA=F)

swissFit$parameters$summary


plot(swissFit$raster[["predict.exp"]]) 
\dontrun{
load(url("http://www.filefactory.com/file/frd1mhownd9/n/CHE_adm0_RData"))
plot(spTransform(gadm, swissRain@proj4string), add=TRUE)
}

}


 	