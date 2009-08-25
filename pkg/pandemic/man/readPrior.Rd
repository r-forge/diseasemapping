\name{readPrior}
\Rdversion{1.1}
\alias{readPrior}
\alias{writePrior}
\title{
Read and write prior distributions from files}
\description{
Hyperparameters are read from files and stored as an R list.
}
\usage{
x=readPrior(file = "priors.txt")
writePrior(x, file = "priors.txt")
}
%- maybe also 'usage' for other objects documented here.
\arguments{
\item{x}{a list of prior distributions}
  \item{file}{
a file name}
}
\seealso{
\code{link[pandemicPriors]{pandemicPriors}}
}
\examples{
x=pandemicPriors()
writePrior(x,"myfile.txt")
x=readPrior("myfile.txt")
}
