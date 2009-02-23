\name{formatCases}
\alias{formatCases}

\title{ Format the case data set }
\description{
 The formatCases funtion cleans and formats the case data set. 
 Changes the other format of age and sex group to two columns: age and sex.
}
\usage{
formatCases(casedata, breaks = NULL, years = NULL)
}
\arguments{
  \item{casedata}{case data set in data.frame format }
  \item{breaks}{ age breaks }
  \item{years}{ number of years}
}
\details{
  ~~ If necessary, more details than the description above ~~
}
\value{
  
}
\references{ }
\author{  }
\note{ 
}
\seealso{  }
\examples{
data(casedata)
cases <- formatCases(casedata, breaks= c(5, 10, 15, 50))

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
