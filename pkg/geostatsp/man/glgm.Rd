\name{glgm-methods}
\docType{methods}
\alias{glgm}
\alias{glgm-methods}
\alias{glgm,ANY,ANY,ANY,ANY-method}
\alias{glgm,formula,Raster,ANY,ANY-method}
\alias{glgm,formula,Spatial,ANY,ANY-method}
\alias{glgm,formula,data.frame,Raster,data.frame-method}


\alias{lgcp}

\title{
Generalized Linear Geostatistical Models
}
\description{
Fits a generalized linear geostatistical model or a log-Gaussian Cox process
using \code{inla} }
\usage{
\S4method{glgm}{ANY,ANY,ANY,ANY}(formula, data,  grid, covariates, buffer=0, shape=1, prior, ...) 
\S4method{glgm}{formula,Raster,ANY,ANY}(formula, data,  grid, covariates, buffer=0, shape=1, prior, ...) 
\S4method{glgm}{formula,Spatial,ANY,ANY}(formula, data,  grid, covariates, buffer=0, shape=1, prior, ...) 
\S4method{glgm}{formula,data.frame,Raster,data.frame}(formula, data,  grid, covariates, buffer=0, shape=1, prior, ...) 
lgcp(formula=NULL, data,  grid, covariates=NULL, border, ...)
}
 
\arguments{
\item{data}{
An object of class \code{ SpatialPointsDataFrame} containing the data.
}
\item{grid}{Either an integer giving the number of cells in the x direction, or a raster object which 
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
 \code{ f(space,..)} term. For \code{glgm} the response variable defaults to the first variable in the \code{data} object, and
 \code{formula} can be an integer or character string specifying the response variable.  For \code{lgcp}, the formula
 should be one-sided.
}
\item{prior}{list with elements named \code{range}, \code{sd}, \code{sdObs}.  See Details.}
\item{shape}{Shape parameter for the Matern correlation function, must be 1 or 2.}
\item{buffer}{Extra space padded around the data bounding box to reduce edge effects.}
\item{border}{boundary of the region on which an LGCP is defined, passed to \code{\link[raster]{mask}}}
\item{...}{Additional options passed to  \code{\link[INLA]{inla}} }
}
\details{
This function performs Bayesian inference for generalized linear geostatistical models with INLA.  The Markov random field
approximation on a regular lattice is used for the spatial random effect.  The range parameter is the distance at which 
the correlation is 0.13, or
\deqn{cov[U(s+h), U(s)] = (2^{1-\nu}/Gamma(\nu))  d^\nu besselK(d, \nu) } 
\deqn{d= |h| \sqrt{8 \nu}/range}
where \eqn{\nu} is the shape parameter. The range parameter produced by \code{glgm} multiplies the range parameter from \code{INLA} by the cell size.

Elements of \code{prior} can be named \code{range}, \code{sd}, or \code{sdObs}.  Elements can consist of: 
\itemize{
\item a single value giving the prior median for penalized complexity priors (exponential on the sd or 1/range).
\item a vector \code{c(u=a, alpha=b)} giving an quantile and probability for pc priors.  For standard deviations alpha is an upper quantile, for the range parameter b = pr(1/range > 1/a).
\item a vector \code{c(lower=a, upper=b)} giving a 0.025 and 0.975 quantiles for the sd or range.
\item a list of the form \code{list(prior='loggamma', param=c(1,2))} passed directly to inla.
\item a two-column matrix of prior densities for the sd or range.
}


}

\value{
A list with two components named  \code{inla}, \code{raster}, and \code{parameters}.  \code{inla} contains the results of the call to the
\code{inla} function.  \code{raster} is a raster stack with the following layers:
\item{random.}{mean, sd, X0.0??quant: Posterior mean, standard deviation, and quantiles of the random effect}
\item{predict.}{mean, sd, X0.0??quant: same for linear predictors, on the link scale}
\item{predict.exp}{posterior mean of the exponential of the linear predictor}
\item{predict.invlogit}{Only supplied if a binomial response variable was used.}

\code{parameters} contains a list with elements:
\item{summary}{a table with parameter estimates and posterior quantiles}
\item{range, sd}{prior and posterior distributions of range and standard deviations}
}

\seealso{
\url{http://r-inla.org}
}
\examples{
# use the 'safe' version of INLA on linux systems
if(Sys.info()['sysname'] =='Linux' &
  requireNamespace("INLA")) {   
  INLA::inla.setOption(inla.call = system.file(paste("bin/linux/",          
          ifelse(.Machine$sizeof.pointer == 4, "32", "64"),
          'bit/inla.static', sep=''), package="INLA")) }

# geostatistical model for the swiss rainfall data
require("geostatsp")
data("swissRain")
swissRain$lograin = log(swissRain$rain)
swissFit =  glgm(formula="lograin", data=swissRain, grid=30, 
	covariates=swissAltitude, family="gaussian", buffer=20000,
	prior = list(sd=1, range=100*1000, sdObs = 2),
	control.inla = list(strategy='gaussian')
	)

if(!is.null(swissFit$parameters) ) {
	
	swissExc = excProb(swissFit, threshold=log(25))

	swissExcRE = excProb(swissFit$inla$marginals.random$space, 
		log(1.5),template=swissFit$raster)

	swissFit$parameters$summary

	plot(swissFit$parameters$range$prior,type="l",
		ylim=c(0,max(swissFit$parameters$range$posterior[,"y"])),
		xlim=c(0, 500000))
	abline(v=swissFit$parameters$range$userPriorCI,col="yellow")
	abline(v=swissFit$parameters$range$priorCI,col="orange")
	lines(swissFit$parameters$range$posterior, col='blue')


}


if(interactive()  | Sys.info()['user'] =='patrick') {
	plot(swissFit$raster[["predict.exp"]]) 

	mycol = c("green","yellow","orange","red")
	mybreaks = c(0, 0.2, 0.8, 0.95, 1)
	plot(swissBorder)
	plot(swissExc, breaks=mybreaks, col=mycol,add=TRUE,legend=FALSE)
	plot(swissBorder, add=TRUE)
	legend("topleft",legend=mybreaks, fill=c(NA,mycol))


	plot(swissBorder)
	plot(swissExcRE, breaks=mybreaks, col=mycol,add=TRUE,legend=FALSE)
	plot(swissBorder, add=TRUE)
	legend("topleft",legend=mybreaks, fill=c(NA,mycol))
}

		

\dontrun{
load(url("http://www.filefactory.com/file/frd1mhownd9/n/CHE_adm0_RData"))
thenames = GNcities(bbox(gadm),max=12)
swissTiles = openmap(bbox(gadm),zoom=8,type="nps")

par(mar=c(0,0,0,0))
plot(gadm)
plot(swissTiles, add=TRUE)
library('RColorBrewer')
mycol=rev(brewer.pal(4,"RdYlGn"))
plot(mask(
		projectRaster(swissExc, crs=proj4string(gadm)),
		gadm), 
	breaks = c(0, 0.2, 0.8, 0.95, 1.00001), 
	col=mycol, alpha=0.5,add=TRUE)	
plot(gadm, add=TRUE, lwd=2, border='blue')

points(thenames,cex=0.5)
text(thenames, labels=thenames$name,pos=3,
  vfont=c("gothic german","plain"),cex=1.5)



}

# a log-Gaussian Cox process example

if(interactive()  | Sys.info()['user'] =='patrick') {
myPoints = SpatialPoints(cbind(rbeta(100,2,2), rbeta(100,3,4)))
myPoints@bbox = cbind(c(0,0), c(1,1))

mycov = raster(matrix(rbinom(100, 1, 0.5), 10, 10), 0, 1, 0, 1)
names(mycov)="x1"


res = lgcp(data=myPoints, grid=20, covariates=mycov,
	formula=~factor(x1),
	priorCI=list(sd=c(0.9, 1.1), range=c(0.4, 0.41))
)
plot(res$raster[["predict.exp"]])
plot(myPoints,add=TRUE,col="#0000FF30",cex=0.5)

}

}


 	