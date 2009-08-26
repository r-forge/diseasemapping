\name{mcmcPandemic}
\Rdversion{1.1}
\alias{mcmcPandemic}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
MCMC simulation for Pandemic
}
\description{
\code{mcmcPandemic} is a function for generating the parameters for different disease stages using the MCMC simulation.
}
\usage{
mcmcPandemic(xdata, params, prior, sigma, runs, thin=1)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{xdata}{the data set simulated using the \code{simEpidemic} function, which contains the informations of number of days in different disease stages for each patients.}
  \item{params}{the parameters simulated using the \code{pandemicParams} function, which contains the esitmated days be in different disease stages.}
  \item{prior}{the priors simualted using the \code{pandemicPriors} function.}
  \item{sigma}{either a scalar or a list.  If a list, it's in the similar format as \code{params}, with each element containing the  standard deviation for random walk Metropolis algorithm for the corresponding parameter. If a scalar, the standard deviation is the starting value given in \code{params} times \code{sigma}.}
  \item{runs}{number of iterations to retain.}
  \item{thin}{number of iterations performed for every iteration retained.  Defaults to 1 for no thinning.}
}
\details{

}
\value{

}
\references{

}
\author{
Patrick Brown
}
\note{

}



\seealso{

}
\examples{

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{}
\keyword{}
