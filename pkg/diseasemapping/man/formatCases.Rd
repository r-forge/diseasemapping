\name{formatCases}
\alias{formatCases}
\title{format the case data set}
\description{
  The formatCases funtion cleans and formats the case data set. Changes the other format of age and sex group to two columns: age and sex.
}
\usage{
formatCases(casedata)
}

\arguments{
  \item{casedata}{case data set in data frame format}
}
\details{
Usually there are two format of age and sex group; 1. F5\_9 or FEMALE5\_9 means female in the age group 5 to 9 years old; 2.
101 or 217 means male in the age group 5 to 9 and female with age 85 plus respectively.  
}
 
\value{
  A dataframe object
}

\author{Patrick Brown}


\examples{
library(diseasemapping)
data(casedata)
casedataclean <- formatCases(casedata)
head(casedataclean)
}

