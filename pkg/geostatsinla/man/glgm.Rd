\name{glgm}
\alias{glgm}
\alias{lgcp}
\alias{geostatsinla-package}

\title{
Generalized Linear Geostatistical Models
}
\description{
Fits a generalized linear geostatistical model or a log-Gaussian Cox process
using   \code{\link[INLA]{inla}} }
\usage{
glgm(data,  cells, covariates=NULL, formula=NULL, priorCI=NULL, 
rough=1, buffer = 0, mesh=FALSE,...)
lgcp(data,  cells, covariates=NULL, formula=NULL, priorCI=NULL, 
rough=1, buffer = 0, mesh=FALSE,...)
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
The spatial random effect should not be supplied but the default 
can be overridden with a 
 \code{ f(space,..)} term. For \code{glgm} the response variable defaults to the first variable in the \code{data} object, and
 \code{formula} can be an integer or character string specifying the response variable.  For \code{lgcp}, the formula
 should be one-sided.
}
\item{priorCI}{list with named elements of 0.025 and 0.975 quantiles of prior distributions, or a single 
vector giving the prior quantiles for the range parameter.  List elements can be named: \code{range} for 
the range parameter (not the scale); \code{sd} for the standard deviation (not the variance or precision); 
\code{intercept} for the intercept; and \code{betas} for the regression coefficients.}
\item{rough}{roughness parameter for the Matern correlation function, must be 1 or 2.}
\item{buffer}{Extra space padded around the data bounding box to reduce edge effects.}
\item{mesh}{Currently unimplemented options for using a mesh rather than a grid for the Markov random field approximation.}
\item{...}{Additional options passed to  \code{\link[INLA]{inla}}}
}
\details{
This function performs Bayesian inference for generalized linear geostatistical models with INLA.  The Markov random field
approximation on a regular lattice is used for the spatial random effect.  The range parameter is the distance at which 
the correlation is 0.13, or
\deqn{cov[U(s+h), U(s)] = 
		(2^(1-rough)/gamma(rough) )*d^rough*besselK(d, rough) } 
\deqn{d= |h|*sqrt(8*rough)/range}
The range parameter produced by \code{glgm} multiplies the range parameter from \code{INLA} by the cell size.
}
\value{
A list with two components named  \code{inla}, \code{raster}, and \code{parameters}.  \code{inla} contains the results of the call to the
\code{\link[INLA]{inla}} function.  \code{raster} is a raster stack with the following layers:
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
\code{\link[INLA]{inla}}
}
\examples{
# geostatistical model for the swiss rainfall data
require(geostatsp)
data(swissRain)
swissRain$lograin = log(swissRain$rain)
swissFit =  glgm(swissRain, cells=30, formula="lograin",
	covariates=swissAltitude, family="gaussian", buffer=20000,
	priorCI=list(sd=c(0.2, 2), range=c(50000,500000)), 
	control.mode=list(theta=c(1.9,0.15,2.6),restart=TRUE),
	control.family=list(hyper=list(prec=list(prior="loggamma", param=c(.1, .1))))
	)

swissFit$parameters$summary

swissExc = excProb(swissFit$inla$marginals.random$space, 0, swissFit$raster)
plot(swissExc, breaks = c(0, 0.2, 0.8, 0.95, 1.00001), 
	col=c('green','yellow','orange','red'))	
plot(swissBorder, add=TRUE)		
		

\dontrun{
load(url("http://www.filefactory.com/file/frd1mhownd9/n/CHE_adm0_RData"))
thenames = GNcities(bbox(gadm),max=12)
swissTiles = openmap(bbox(gadm),zoom=8,type="nps")

par(mar=c(0,0,0,0))
plot(gadm)
plot(swissTiles, add=TRUE)
library(RColorBrewer)
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

plot(swissFit$parameters$range$prior,type="l",
	ylim=c(0,max(swissFit$parameters$range$posterior[,"y"])),
	xlim=c(0, swissFit$parameters$range$userPriorCI[2]*1.1))
abline(v=swissFit$parameters$range$userPriorCI,col="yellow")
abline(v=swissFit$parameters$range$priorCI,col="orange")
lines(swissFit$parameters$range$posterior, col='blue')

plot(swissFit$raster[["predict.exp"]]) 
# a log-Gaussian Cox process example

myPoints = SpatialPoints(cbind(rbeta(100,2,2), rbeta(100,3,4)))
myPoints@bbox = cbind(c(0,0), c(1,1))

mycov = raster(matrix(rbinom(100, 1, 0.5), 10, 10), 0, 1, 0, 1)
names(mycov)="x1"

res = lgcp(data=myPoints, cells=20, covariates=mycov,
	formula=~factor(x1),
	priorCI=list(sd=c(0.9, 1.1), range=c(0.4, 0.41))
)
plot(res$raster[["predict.exp"]])
plot(myPoints,add=TRUE,col="#0000FF30",cex=0.5)

# a model with little data, posterior should be same as prior

data2 = SpatialPointsDataFrame(cbind(c(1,0), c(0,1)), data=data.frame(y=c(NA,NA), offset=c(-100,-200)))

res = glgm(data2, cells=20, formula=y~1, priorCI = list(sd=c(1,2), range=c(0.3, 2)),
family="poisson",buffer=2)

priorPrec = res$par$sd$params
priorRange = res$par$range$params
 par(mfrow=c(2,1))
# sd
plot(res$parameters$sd$prior,type='l', col='blue')
lines(res$parameters$sd$post,col='red')


# range
plot(res$parameters$range$prior,type='l', col='blue')
lines(res$parameters$range$post,col='red')
legend("topright", col=c("blue","red"),lty=1,legend=c("prior","post'r"))

}


 	