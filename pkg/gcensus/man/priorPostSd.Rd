\name{priorPostSd}
\alias{priorPostSd}
\alias{precToSd}
\title{precisions to standard deviations}
\usage{
priorPostSd(res, param = 1:length(res$all.hyper$random), minSd = 0.001)
precToSd(densmat)
}
\arguments{
\item{res}{an \code{\link[INLA]{inla}} result}

\item{param}{vector of parameters to transform}

\item{minSd}{minimum value of standard deviation to consider}
\item{densmat}{matrix with two columns with precisions and densities}
}
\description{
Transforms prior and posterior distributions of precision parameters to standard deviations
}

