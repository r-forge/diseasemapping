\name{area}
\alias{area}

\title{calculate the area of the region}
\description{
  Take a single spatial polygon object and calculate the area of each region, return a vector of areas.
  
}
\usage{
area(sp)
}

\arguments{
  \item{sp}{name of a spatial polygon data frame.}
}


\author{ Patrick Brown }


\examples{
data(popdata)
area(popdata)
}

