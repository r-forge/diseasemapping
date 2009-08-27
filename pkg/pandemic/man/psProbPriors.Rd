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
psPrior( taub1=2,
 taub2=0.05,
 priorMean=0.01,
 upper95 = 0.02)
}
\arguments{
  \item{fatality}{
Prior for fatality probabilities
}
  \item{hosp}{
Prior for hospitalization probabilities (conditional on not being fatal) 
}
  \item{taub1}{
Shape parameter for gamma prior on precision of spline coefficients}
  \item{taub2}{
Range parameter for gamma prior on precision of spline coefficients}
  \item{priorMean}{
Prior mean for the baseline probability}
  \item{upper95}{
Prior upper 95 pct boundary for baseline probability}
}
\value{
 \code{psPrior} returns a list with elements taub1, taub2, beta0 (prior mean on logit scale), and Sbeta0 (prior variance of logit mean).
 
\code{psProbPriors} returns a list with elements \code{fatality} and \code{hosp}, each of which are returned by \code{psPrior}.
}
\seealso{
\code{\link[DPpackage:PSgam]{PSgam}}
}
\examples{
priors= pandemicPriors(probs=psProbPriors(hosp=psPrior(tau1=1)))
priors$fatality
}
