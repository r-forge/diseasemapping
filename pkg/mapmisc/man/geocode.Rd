\name{geocode}
\alias{geocode}
\title{Georeferencing with Google}
\description{
Uses the dismo package to geocode with Google
  }
\usage{
	geocode(...)
	}

\arguments{
	\item{...}{arguments passed to \code{\link[dismo]{geocode}}}
}

\seealso{\code{\link[dismo]{geocode}}}


\value{A \code{SpatialPointsDataFrame} in long-lat projection.}


\examples{


\dontrun{

if (requireNamespace("dismo", quietly = TRUE)) { 

cities=geocode('Ulan batar')
mytiles = openmap(cities, buffer=800*1000)

map.new(mytiles)
plot(mytiles, add=TRUE)
points(cities, col='red')
text(cities, labels=cities$originalPlace, col='red',pos=4)

}
}


}