\name{getBreaks}
\alias{getBreaks}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ ~~function to do ... ~~ }
\description{
  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
getBreaks(colNames, breaks = NULL, mustAggregate = TRUE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{colNames}{ ~~Describe \code{colNames} here~~ }
  \item{breaks}{ ~~Describe \code{breaks} here~~ }
  \item{mustAggregate}{ ~~Describe \code{mustAggregate} here~~ }
}
\details{
  ~~ If necessary, more details than the description above ~~
}
\value{
  ~Describe the value returned
  If it is a LIST, use
  \item{comp1 }{Description of 'comp1'}
  \item{comp2 }{Description of 'comp2'}
  ...
}
\references{ ~put references to the literature/web site here ~ }
\author{ ~~who you are~~ }
\note{ ~~further notes~~ 

 ~Make other sections like Warning with \section{Warning }{....} ~
}
\seealso{ ~~objects to See Also as \code{\link{help}}, ~~~ }
\examples{
##---- Should be DIRECTLY executable !! ----
##-- ==>  Define data, use random,
##--	or do  help(data=index)  for the standard data sets.

## The function is currently defined as
function (colNames, breaks=NULL, mustAggregate = TRUE) 
{
    popColumns <- grep("^(m|f|male|female)[[:digit:]]+(_|-|plus|\\+)[[:digit:]]*$", 
        colNames, value = T, ignore.case = T)
    ageGroups <- gsub("^(m|f|male|female)", "", popColumns, ignore.case = TRUE)
    ageGroups <- gsub("(\\+|plus)", "_Inf", ageGroups, ignore.case = TRUE)
    ageLower <- as.numeric(gsub("(_|-)([[:digit:]]+|Inf)$", "", 
        ageGroups))
    ageUpper <- as.numeric(gsub("^[[:digit:]]+(_|-)", "", ageGroups))
    currentbreaks <- c(sort(unique(ageLower)), max(ageUpper))
    sex <- substr(popColumns, 1, 1)

    result = list(breaks = currentbreaks, age = ageLower, 
        sex = sex, oldNames = popColumns, newNames = paste(sex, 
            ageLower, sep = "."), mustAggregate = mustAggregate)

if(length(breaks)> 0 ) {
 # check that the breaks are valid

# make sure that the right end of the last bin is Inf
	if(breaks[length(breaks)] != Inf) {
	  if(breaks[length(breaks)] \%in\% currentbreaks) {
		 breaks = c(breaks, Inf)
	  } else {
		breaks[length(breaks)] = Inf
	  }
	}

# if the breaks don't include some younger age groups, keep those age groups as-is
breaks = c(currentbreaks[currentbreaks < min(breaks)], breaks)

# check that population breaks are nested within breaks
if (all(breaks \%in\% currentbreaks)) {

# if the breaks are different from the population breaks
if(any(breaks != currentbreaks)) {
 result$breaks = breaks
 result$mustAggregate = F
  }

  } else {
	warning("population data doesn't appear to be nested within breaks, ignoring breaks")
  }
  }
result
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
