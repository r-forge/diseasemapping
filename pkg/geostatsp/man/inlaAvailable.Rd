\name{inlaAvailable}
\alias{inlaAvailable}
\alias{inlaSetThreads}
\title{Check whether INLA is usable}
\description{
Tests whether the INLA package is installed and its binaries work, and
optionally sets INLA thread options.
}
\usage{
inlaAvailable()
inlaSetThreads(num.threads = 2L)
}
\arguments{
\item{num.threads}{Number of threads passed to
\code{INLA::inla.setOption}.}
}
\details{
\code{requireNamespace("INLA")} can succeed even when the INLA program
binaries are missing or broken.  \code{inlaAvailable} additionally
probes \code{INLA::inla.getOption} so examples and tests can skip safely.
\code{inlaSetThreads} is a no-op when INLA is unavailable.
}
\value{
\code{inlaAvailable} returns a logical.  \code{inlaSetThreads} returns
\code{TRUE} if options were set, \code{FALSE} otherwise, invisibly.
}
\examples{
inlaAvailable()
inlaSetThreads(2)
}
