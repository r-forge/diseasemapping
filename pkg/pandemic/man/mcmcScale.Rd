\name{mcmcScale}
\Rdversion{1.1}
\alias{mcmcScale}

\title{
Create standard deviations for metropolis hastings updating from starting values
}
\description{
Defines standard deviations to be the starting values times the supplied scaling factor, subject to maximums and minimums.
}
\usage{
mcmcScale(params, sigma, minScale = 0.01, maxScale = 0.2)
}
\arguments{
  \item{params}{
starting values}
  \item{sigma}{
the scaling factor}
  \item{minScale}{
minimum possible standard deviation}
  \item{maxScale}{
maximum standard deviation}
}
\value{
A list in the same format as the parameters, containing standard deviations instead of parameter values}
