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
\item{extent}{an Extent object, or any object from which an Extent can be obtained.}
\item{progress}{see \code{\link[dismo]{geocode}}}
\item{...}{additional arguments passed to \code{\link[dismo]{geocode}}}
}

\seealso{\code{\link[dismo]{geocode}}}


\value{A \code{SpatialPointsDataFrame} with coordinates in the 
projection of \code{extent} if possible, or long-lat
otherwise.}

\details{If the option  
\code{options()$mapmiscCachePath} is set, it will be used
to specify the folder to save downloaded 
data.}

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