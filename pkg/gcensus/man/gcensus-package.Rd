\name{gcensus}
\alias{gcensus}
\alias{gcensus-package}
\title{Retrieve and combine us boundary, population and income data}
\usage{
gcensus(country = "USA", year = 2010, level = 1, subset = NULL,
  database = c(dbname = "gcensus", user = "gcensus", host = "localhost", port
  = 5433), income = TRUE, crs = "EPSG:4326", simplifyTolerance = 0)
}
\arguments{
\item{country}{string specifying a country}

\item{year}{census year}

\item{level}{see below}

\item{subset}{passed as \code{where} to postgreql}

\item{database}{details of database to retrieve data from}

\item{income}{logical, retrieve income data}

\item{crs}{coordinate reference system for the output}

\item{simplifyTolerance}{if non-zero, simplify polygons with this tolerance}
}
\description{
Retrieve and combine us boundary, population and income data
}
\details{
level	abbreviation	name
0	USA	United States of America
1	STATE	State
1.1	PUMA	Public Use Microdata Areas
1.2	ZCTA	Zip Code Tabulation Areas
2	COUNTY	County
2.1	COSUB	County Subdivision
3	TRACT	Tract
4	BLKGRP	Block Group

Canadian Census:

level	abbreviation	name
0	USA	United States of America
1	STATE	State
1.1	PUMA	Public Use Microdata Areas
1.2	ZCTA	Zip Code Tabulation Areas
2	COUNTY	County
2.1	COSUB	County Subdivision
3	TRACT	Tract
4	BLKGRP	Block Group

}

\examples{
\dontrun{
usa = gcensus(country = "USA", year = 2010, level = 1, subset = NULL,
  database = c(dbname = "gcensus", user = "gcensus", 
  host = "localhost", port = 5433), 
  income = TRUE, crs = "EPSG:4326", simplifyTolerance = 0)
}
}