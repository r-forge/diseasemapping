\name{plotAgeProbs}
\Rdversion{1.1}
\alias{plotAgeProbs}

\title{
Plot posterior distribution of age dependent probabilities
}
\description{
Plots posterior mean, quantiles and (if provided) prior quantiles.
}
\usage{
plotAgeProbs(postSample, type = "fatality", quantiles = c(0.025, 0.975), prior = NULL)
}
\arguments{
  \item{postSample}{
sample from the posterior distribution, produced by \code{mcmcPandemic}
}
  \item{type}{
which probability to plot, fatality or hospital.
}
  \item{quantiles}{
Quantiles of the posterior distribution to plot
}
  \item{prior}{
the prior distribution used in the analysis.
}
}

