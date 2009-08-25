\name{psProbPriors}
\Rdversion{1.1}
\alias{psProbPriors}
\alias{psPrior}
\title{
Priors for PSgam
}
\description{
Create prior distributions for age-dependent probabilities in pandemic models
}
\usage{
psProbPriors(fatality = psPrior(), hosp = psPrior())
psPrior(taub1 = 0.01, taub2 = 0.01, beta0 = rep(0, 1), Sbeta0 = diag(100, 1), tau1 = 0.01, tau2 = 0.01)
}
\arguments{
  \item{fatality}{
Prior for fatality probabilities
}
  \item{hosp}{
Prior for hospitalization probabilities (conditional on not being fatal) 
}
  \item{taub1}{
PSgam parameter}
  \item{taub2}{
PSgam parameter}
  \item{beta0}{
PSgam parameter}
  \item{Sbeta0}{
PSgam parameter}
  \item{tau1}{
PSgam parameter}
  \item{tau2}{
PSgam parameter}
}
\value{
 \code{psPrior} returns a list with all the function arguments as elements. 
\code{psProbPriors} returns a list with elements \code{fatality} and \code{hosp}, each of which are returned by \code{psPrior}.
}
\seealso{
\code{\link{DPpackage:PSgam}}}
\examples{
priors= pandemicPriors(probs=psProbPriors(hosp=psPrior(tau1=1)))
priors$fatality
}