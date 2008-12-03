\name{formatPopulation}
\alias{formatPopulation}
\alias{formatPopulation.data.frame}
\alias{formatPopulation.list}
\alias{formatPopulation.PostgreSQLConnection}
\alias{formatPopulation.SpatialPolygonsDataFrame}

\title{ Format and clean the population data set}
\description{
The formatCases funtion cleans and formats the population data set. Changes the spatial polygon shape file to data frame 

}
\usage{
formatPopulation(popdata, aggregate.by, ...)
formatPopulation.data.frame(popdata, aggregate.by = NULL, ...)
formatPopulation.list(popdata, aggregate.by = NULL, cyears = as.integer(names(popdata)),  times = "YEAR", ...)
formatPopulation.PostgreSQLConnection(popdata, aggregate.by = NULL, subset = NULL, ...)
formatPopulation.SpatialPolygonsDataFrame(popdata, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{popdata}{  population data set, either a list, data frame, sp object, or database connection.}
  \item{aggregate.by}{the term that used to aggregate the number of people. i.e year }
  \item{cyears}{vector of census year corresponding to the elements of the popdata list}
  \item{subset}{argument for selecting a subset of the data}

  \item{times}{ the variable name which stores time }
 \item{\dots}{ additinal infomation }
 }
\details{
 This function reshape the data frame in the function in "long" format which could aggregate the number of population according to the age and set group.
}
\value{

}
\references{  }
\author{ }
\note{ 
}
\seealso{\code{\link[pkg:diseasemapping]{formatPopulation.data.frame}} from the \code{diseasemapping} pagkage, 
         \code{\link[pkg:diseasemapping]{formatPopulation.list}}, 
         \code{\link[pkg:diseasemapping]{formatPopulation.PostgreSQLConnection}},
         \code{\link[pkg:diseasemapping]{formatPopulation.SpatialPolygonsDataFrame}}. 
 }
\examples{
# if population data set is a shape file.
data(popdata)
head(popdata@data)
popdatatry <- formatPopulation(popdata)
head(popdatatry)

# if population data set is a data frame
formatPopulation.data.frame(popdata@data)

# if population data set is a list
#formatPopulation.list(popdata@data)

# if population data set is PostgreSQLConnection
#formatPopulation.PostgreSQLConnection(popdata)

# if population data set is SpatialPolygonsDataFrame
formatPopulation.SpatialPolygonsDataFrame(popdata)

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ }
\keyword{ }% __ONLY ONE__ keyword per line
