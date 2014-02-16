\name{loglikGmrf}
\alias{loglikGmrf}
\alias{loglikGmrfNugget}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Likelihood of a Gaussian Markov Random Field model
}
\description{
Calculates the likelihood of a Matern-approximating GMRF on 
a regular grid.
}
\usage{
loglikGmrfNugget(ar, Yvec, Xmat, NN, maternShape = 1)
loglikGmrfNugget(ar, propNugget = 5, Yvec, Xmat, NN, maternShape = 1)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{ar}{
An autoregressive parameter between 0 and 1, see Details.
}
  \item{propNugget}{
Nugget effect variance divided by the conditional variance of the spatial process.
}
  \item{Yvec}{
Vector of observations
}
  \item{Xmat}{
Design matrix of covariates (including intercept and with indicator variables)
}
  \item{NN}{
Nearest neighbour matrix as produced by  \code{\link[geostatsp]{NNmat}}
}
  \item{maternShape}{
order of the Matern correlation function.
}
}
\details{
The weights of the GMRF sum to one and ar multiplied by the \code{ar} argument.
The resulting GMRF approximates a Matern correlation 
function with range parameter of 
\code{sqrt(ar/(1-ar))*sqrt(2*maternShape)} grid cells, using the parametrization
in  \code{\link[geostatsp]{matern}}.  
}
\value{
A vector of values including the likelihood, conditional variance, marginal variance,
and range. 
}

\author{
Patrick Brown
}


\seealso{
 \code{\link[geostatsp]{NNmat}},  \code{\link[geostatsp]{matern}}
}
\examples{

data('swissRainR')

swissRainDf = as.data.frame(swissRainR)
Xmat = 	model.matrix( ~ alt,swissRainDf)
swissShape = 1

loglikGmrf(0.8, swissRainDf[,'prec5'],
	Xmat, swissNN, swissShape)

loglikGmrfNugget(ar=0.8, propNugget=0.5,
	Yvec=swissRainDf[,'prec5'],
	Xmat=Xmat, NN=swissNN, maternShape=swissShape)


}


