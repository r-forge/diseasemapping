\name{priorPost}
\alias{priorPostSd}
\alias{priorPost}
\alias{precToSd}
\title{precisions to standard deviations}
\usage{
priorPost(object)
priorPostSd(object, param,
  group = c("random", "family"))
precToSd(densmat)
}
\arguments{
\item{object}{an \code{inla} result}

\item{param}{vector of parameters to transform}

\item{group}{random effect parameters or 'family' parameters}

\item{densmat}{matrix with columns 'x' and 'y'}
}
\description{
Transforms prior and posterior distributions of precision parameters to standard deviations
}
