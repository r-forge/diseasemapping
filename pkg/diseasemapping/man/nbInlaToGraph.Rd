\name{nbToInlaGraph}
\alias{nbToInlaGraph}


\title{
Write a graph file for INLA
}
\description{
Writes a graph file from an adjacency matrix suitable for use with INLA.
}

\usage{
nbToInlaGraph(adjMat, graphFile="graph.dat")
}

\arguments{
\item{adjMat}{An object of class \code{nb} containing the adjacency matrix.}
\item{graphFile}{name of file to save adjacency matrix to.}
}

\value{
	A vector of names and indices
}

\author{
Patrick Brown
}


\seealso{
\code{\link[spdep]{poly2nb}},
\code{\link[spdep]{nb2INLA}} 
}

\examples{

data('kentucky')

if( require("spdep", quietly=TRUE)) {

kentuckyAdj = poly2nb(
	kentucky, row.names=kentucky$County,
	queen=FALSE
)

attributes(kentuckyAdj)$region.id[1:5]
unlist(lapply(kentuckyAdj, length))[1:5]

nFile = tempfile()

nbRes = nbToInlaGraph(kentuckyAdj, nFile)

nbRes[1:5]

readLines(nFile, n=6)

}

}

