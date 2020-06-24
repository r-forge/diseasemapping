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


cities=try(mapmisc::geocode('Ulan batar'), silent=TRUE)
data('worldMap', package='mapmisc')

if(
	!all(class(cities) == 'try-error') & requireNamespace("rgdal", quietly = TRUE)
	) { 
citiesT = spTransform(cities, crs(worldMap))
mapmisc::map.new(citiesT, buffer=5000*1000)
plot(worldMap, add=TRUE)
points(citiesT, col='red')
text(citiesT, labels=citiesT$name, col='red',pos=4)

}



}