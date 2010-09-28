\name{simEpidemic}
\Rdversion{1.1}
\alias{simEpidemic}
\title{
Data Simulation
}
\description{
Using the parameters generated from the \code{pandemicParams} function, then simulate the number of patients for different stages and the case levels.
}
\usage{
simEpidemic(params, delta = 5, days = 20, probOnsetMissing = 0.7, randomInfections = TRUE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{params}{the parameters generated from the \code{pandemicParams} function.}
  \item{delta}{mean number of infections per day}
  \item{days}{number of study days}
  \item{probOnsetMissing}{the probability of missing information for onset cases.}
  \item{randomInfections}{a logical indicating whether the patients randomly get infection. Default is TRUE. }
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
someParams = pandemicParams()
somedata = simEpidemic(someParams, delta=10)

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{  }
\keyword{  }% __ONLY ONE__ keyword per line
                                            