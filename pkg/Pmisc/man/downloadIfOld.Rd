\name{downloadIfOld}
\alias{downloadIfOld}
\title{Download a file}
\usage{
downloadIfOld(url, file = basename(url), age = "5 days", 
	verbose = FALSE,  exdir=tempdir())
}
\arguments{
\item{url}{a string specifying the remote location}

\item{file}{local file name}

\item{age}{maximum age of the local file}

\item{verbose}{print additional information}
\item{exdir}{directory to put contents of any zip files}
}
\value{
A character string of files downloaded or unzipped.
}

\description{
Checks if a downloaded file is old and re-downloads is necessary.  
Any zip files are unzipped.
}

\examples{



\dontrun{

theFiles = Pmisc::downloadIfOld(
'https://cran.r-project.org/src/base/Historic/Windows/mva.zip'
)
theFiles
}

}