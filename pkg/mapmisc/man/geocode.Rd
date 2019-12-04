\name{geocode}
\alias{geocode}
\title{Georeferencing with Google}
\description{
Uses the dismo package to geocode with Google
  }
\usage{
	geocode(x, extent,
	  lang = gsub("(_|[:]).*", "", Sys.getenv('LANGUAGE')))
	}

\arguments{
\item{x}{Vector of character strings to search for}
\item{extent}{Currently unused.  an Extent object, or any object from which an Extent can be obtained.}
\item{lang}{Language for place names in result.}
}



\value{A \code{SpatialPointsDataFrame} with coordinates in the 
projection of \code{extent} if possible, or long-lat
otherwise.}

\details{If the option  
\code{getOption('mapmiscCachePath')} is set, it will be used
to specify the folder to save downloaded 
data. \code{getOption('mapmiscVerbose')} for printing progress.

Data are retreived from Openstreetmap.org, see \url{https://wiki.openstreetmap.org/wiki/Nominatim}.
}

\examples{



if(requireNamespace("rgdal", quietly = TRUE)) { 

cities=try(geocode('Ulan batar'), silent=TRUE)
if(!all(class(cities) == 'try-error')) {
mytiles = openmap(cities, buffer=800*1000, maxTiles=1)

map.new(mytiles)
plot(mytiles, add=TRUE)
points(cities, col='red')
text(cities, labels=cities$name, col='red',pos=4)
}

}


}