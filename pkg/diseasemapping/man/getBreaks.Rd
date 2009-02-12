\name{getBreaks}
\alias{getBreaks}

\title{Calculate the age groups  }
\description{
A function seperate the age group according to the user's wish if it is different from the default age group from the population data set.
}
\usage{
getBreaks(colNames, breaks = NULL)
}

\arguments{
  \item{colNames}{ the names of the columns contains the age info. usually use the column name of the population data set. }
  \item{breaks}{ vector of the age cuts.}
}
\details{
usually, the default age cut from population data set is 5 years group. But sometimes, the user may want to do a 10 years
age breaks. For the older population group which contains less infomation, it might be better to do a 20 year age break.
}
\value{

}
\references{  }
\author{ Lutong Zhou}
\note{ The breaks can only be the a multiple of the default age group. i.e if the default age group is 5 years, 
this function allow the age breaks change to 10 years or 15 years and so on. But 7 years or 9 years will now be allowed. 
}
\seealso{  }
\examples{
data(popdata)
oldbreak <- getBreaks(names(popdata), breaks=NULL)
newbreak <- getBreaks(names(popdata), breaks=c(10, 20, 50, 70))

}

\keyword{ }
\keyword{ }% __ONLY ONE__ keyword per line
