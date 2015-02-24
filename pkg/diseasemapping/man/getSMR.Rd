\name{getSMR}
\alias{getSMR}
\title{Calculate the standardized mortality/morbidity ratios   }
\description{
The getSMR function calculates the rate of observe value over expected value. 
It will also merge back the observed value, expected value and the ratio back to the population data set. 
}
\usage{
getSMR(popdata, model, casedata, regionCode , regionCodeCases , 
area = FALSE, area.scale = 1,  ...)
}
\arguments{
  \item{popdata}{ the name of population data set  }
  \item{model}{rates, either fitted model (usually a \code{glm} object), or a vector of rates. }
  \item{casedata}{ the name of case data set }
  \item{regionCode}{ the name of district area column in population data set}
  \item{regionCodeCases}{the name of district area column in case data set }
  \item{area}{ if \code{TRUE}, calculate the expected number of cases per unit area in each region }
  \item{area.scale}{ scale the unit of area. e.g $10^6$: if your spatial coordinates are metres 
  and you want intensity in cases per km2 }
 \item{\dots}{ additional arguments if \code{popdata} is a list of populations from different years.
  these could include \code{years} (a vector with the year of each dataset), \code{year.range} (two dimensional
  vector with first and last year),  
  and \code{personYears} (logical, if \code{TRUE} expected counts are per person year).}

}
\details{
If \code{model} is numeric, it's assumed to be a vector of rates, with the names of the elements corresponding to columns of the population data set.  Names do not need to match exactly (can have M in one set of names, male in another for instance).

Otherwise, \code{model} is passed to the \code{predict} function.  
}
\value{
  Returns a new population data set contains expected number of cases, observed number of cases and SMR. 
  It has the same format as the population data set which put into the function.
}
 
\examples{
data(kentucky)

kentucky2 = getSMR(kentucky, larynxRates, larynx, 
		regionCode="County")

data.frame(kentucky2)[1:10,grep("^F|^M", names(kentucky2), invert=TRUE)]

}
 