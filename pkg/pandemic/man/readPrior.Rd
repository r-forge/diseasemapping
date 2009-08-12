\name{readPrior}
\Rdversion{1.1}
\alias{readPrior}
\aliad{writePrior}
\title{
Read and write prior distributsions from files}
\usage{
x=readPrior(file = "priors.txt")
}
\usage{
writePrior(x, file = "priors.txt")
}
%- maybe also 'usage' for other objects documented here.
\arguments{
\item{x}{a list of prior distributions}
  \item{file}{
a file name}
}
\seealso{
\code{link{pandemicPriors}}
}
\examples{
x=pandemicPriors()
writePrior(x,"myfile.txt")
x=readPriors("myfile.txt")
}
