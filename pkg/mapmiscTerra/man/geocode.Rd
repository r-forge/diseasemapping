\name{geocode}
\alias{geocode}
\title{Georeferencing with OpenStreetMap Nominatim}
\description{
Geocode place names using the OpenStreetMap Nominatim service.
  }
\usage{
	geocode(x, extent,
	  lang = gsub("(_|[:]).*", "", Sys.getenv('LANGUAGE')))
	}

\arguments{
\item{x}{Vector of character strings to search for}
\item{extent}{Currently unused.  An extent, or any object from which an extent can be obtained.}
\item{lang}{Language for place names in result.}
}



\value{A \code{SpatVector} of points with coordinates in the 
projection of \code{extent} if possible, or long-lat
otherwise.}

\details{If the option  
\code{getOption('mapmiscCachePath')} is set, it will be used
to specify the folder to save downloaded 
data. \code{getOption('mapmiscVerbose')} for printing progress.

Data are retrieved from OpenStreetMap.org; see \url{https://wiki.openstreetmap.org/wiki/Nominatim}.
Please respect Nominatim's usage policy when making many requests.
}

\examples{


cities=try(geocode('Ulan batar'), silent=TRUE)
data('worldMap')
worldMap = terra::unwrap(worldMap)

if(!all(class(cities) == 'try-error') ) { 
citiesT = project(cities, crs(worldMap))
oldpar=map.new(citiesT, buffer=5000*1000)
plot(worldMap, add=TRUE)
points(citiesT, col='red')
suppressWarnings(text(citiesT, labels=citiesT$name, col='red',pos=4))
\dontrun{
	# uses unicode symbols
	text(citiesT, labels=citiesT$display_name, col='red',pos=1))
}
par(oldpar)
}



}
