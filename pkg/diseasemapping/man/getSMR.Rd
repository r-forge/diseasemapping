\name{getSMR}
\alias{getSMR}

\title{Calculate the standardized mortality ratios  }
\description{
 The getSMR function calculates the rate of observe value over expected value. It will also merge back the observed value, expected value and the ratio back to the population data set. 
}
\usage{
getSMR(model, popdata, casedata = NULL, regionCode = "CSDUID", regionCodeCases = "CSD2006", cyears = NULL, year.range = NULL, area = FALSE, area.scale = 1)
}

\arguments{
  \item{model}{ the estimated coefficients we calculated from getRates function. }
  \item{popdata}{the population data set name }
  \item{casedata}{the case data set name}
  \item{regionCode}{the name of district area column in population data set }
  \item{regionCodeCases}{ the name of district area column in case data set}
  \item{cyears}{ a vector of census years  }
  \item{year.range}{ study period: a vector of two elements, starting dates and ending dates }
  \item{area}{ if \code{TRUE}, calculate the expected number of cases per area in each region  }
  \item{area.scale}{ control the unit of area. e.g $10^6$: per square kilometers }
}
\details{
  
 If there is no rate calculates from the getRates function, it will leave NA in the merged population data set. 
}
\value{
  Returns a new population data set contains expected number of cases, observed number of cases and SMR. It has the same format as the population data set which put into the function.
}

\author{Patrick Brown }


\examples{
data(casedata)
data(popdata)
model <- getRates(casedata, popdata, ~age*sex, S=c("M","F"))
popdatanew <- getSMR(model, popdata, casedata)
head(popdatanew@data)
}

