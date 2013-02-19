\name{glgm}
\alias{glgm}

\title{
Generalized Linear Geostatistical Models
}
\description{
Fits a generalized linear geostatistical model
using   \code{\link[inla]{inla}} }
\usage{
glgm(formula, data, cells, covariates, prior, ...)
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
will be used for the spatial random effect.}
\item{covariates}{Either a list of rasters or a raster stack containing covariate values used when 
making spatial predictions.  Names of the raster layers correspond to names in the formula.  If
a covariate is missing from the data object it will be extracted from the rasters.}
}
\item{prior}{Prior distribution for the range parameter of the spatial random effect, in units of distance
not numbers of cells.}
}
\details{
This function is essentially a wrapper for   \code{\link[geoR]{krige.conv}}.
}
\value{
If the prediction locations are identical for the fixed effects and random effects, a single raster stack is 
returned with the following layers:
\item{fixed}{Estimated means from the fixed effects portion of the model}
\item{random}{Predicted random effect}
\item{krige.var}{Conditional variance of predicted random effect}
\item{predict}{Prediction of the response, sum of fixed and random effects.  If exp.pred is TRUE, gives predictions
on the exponentated scale, and half of krige.var is added prior to exponentiating}
\item{predict.log}{If exp.pred=TRUE, the prediction of the logged process.}
If the prediction locations are different for fixed and random effects (typically coarser for the random 
effects), a list with two raster stacks is returned.
\item{prediction}{A raster stack as above, though the random effect prediction is resampled to the 
same locations as the fixed effects.}
\item{random}{the predictions and conditional variance of the random effects, on the same
raster as \code{locations}}


}


\seealso{
\code{\link[geoR]{krige.conv}}
}
\examples{
data(swissRain)
swissRain$lograin = log(swissRain$rain)
swissRain$elevation = extract(swissAltitude, swissRain)
swissFit =  likfit(swissRain, lograin ~ elevation, 
	cov.model="matern", kappa = 2, fix.kappa=TRUE,
	ini.cov.pars= c(0.6, 7000))
swissRaster = raster(extent(swissRain@bbox), ncols=61, nrows=91, crs=swissRain@proj4string)	
swissKrige = krige(swissFit, swissRain, swissRaster, swissAltitude)


swissFit =  likfit(swissRain, rain ~ elevation, 
	cov.model="matern", kappa = 2, fix.kappa=TRUE,
	lambda=0.5, fix.lambda=TRUE,
	ini.cov.pars= c(0.6, 7000))
swissKrige = krige(obj.model=swissFit, geodata=swissRain, 
	locations=swissRaster, covariates=swissAltitude)


plot(swissKrige$prediction[["predict"]]) 
\dontrun{
load(url("http://www.filefactory.com/file/frd1mhownd9/n/CHE_adm0_RData"))
plot(spTransform(gadm, swissRain@proj4string), add=T)
}

}

 

