\name{addAgeProbs}
\Rdversion{1.1}
\alias{addAgeProbs}
\title{
Initialize age-dependent probabilities}
\description{
Creates ages and probabilities suitable for pandemic simulation and inference
}
\usage{
addAgeProbs(age = 0:100, prob = rep(0.1, length(age)))
}
\arguments{
  \item{age}{
vector of ages}
  \item{prob}{
vector of probabilities, of the same length as \code{age}}
}
\value{
A data frame with columns age and prob.
}
\examples{
ageProbs = addAgeProbs(age=seq(10, 80, len=8), prob = seq(0.01, 0.1, len=8))
params = pandemicParams(
  probs=NULL,
  ageProbs=list(
    S=addAgeProbs(),
    D=addAgeProbs()
  )
)
}
