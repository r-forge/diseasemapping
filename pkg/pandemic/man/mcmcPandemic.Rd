\name{mcmcPandemic}
\Rdversion{1.1}
\alias{mcmcPandemic}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
MCMC simulation for Pandemic
}
\description{

}
\usage{
mcmcPandemic(xdata, params, prior, sigma, runs)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{xdata}{the data set simulated using the \code{simEpidemic} function, which contains the informations of number of days in different disease stages for each patients.}
  \item{params}{the parameters simulated using the \code{pandemicParams} function, which contains the esitmated days be in different disease stages.}
  \item{prior}{the priors simualted using the \code{pandemicPriors} function.}
  \item{sigma}{the standard deviation for random walk Metropolis algorithm.}
  \item{runs}{number of iterations.}
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
