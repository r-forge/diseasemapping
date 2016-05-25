\name{geocode}
\alias{geocode}
\title{Georeferencing with Google}
\description{
Uses the dismo package to geocode with Google
  }
\usage{
	geocode(x, oneRecord=FALSE, extent=NULL, progress='', ...)
	}

\arguments{
\item{x}{vector of character strings to search for}
\item{oneRecord}{see \code{\link[dismo]{geocode}}}
\item{extent}{}
\item{progress}{see \code{\link[dismo]{geocode}}}
	\item{...}{additional arguments passed to \code{\link[dismo]{geocode}}}
}

\seealso{\code{\link[dismo]{geocode}}}


\value{A \code{SpatialPointsDataFrame} with coordinates in the 
projection of \code{extent} if possible, or long-lat
otherwise.}

\details{If the options \code{options()$mapmisc$verbose} or 
options()$mapmisc$cachePath} are set, they will be used
to control the information printed and folder to save downloaded 
data respectively.}

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