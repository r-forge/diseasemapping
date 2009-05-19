\name{plotPrior}
\Rdversion{1.1}
\alias{plotPrior}

\title{
Plot Priors
}
\description{
Plot the prior distributions to find better priors.
}
\usage{
plotPrior(x, posteriorSample = NULL, file = NULL, quantiles = c(0.025, 0.975), tex = FALSE, transition = "HospDeath")
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{the R object to be plotted, usually it is a list.}
  \item{posteriorSample}{the distribution name of posterior sample.}
  \item{file}{give a name of the file if want to save the priors.}
  \item{quantiles}{confidence interval or range of the paramters. }
  \item{tex}{if \code{TRUE} it will produce a R sweave file for latex.}
  \item{transition}{the name of the stage those prior distributions will be plotted. }
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
\keyword{ }
\keyword{ }% __ONLY ONE__ keyword per line
