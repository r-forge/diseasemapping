\name{downloadIfOld}
\alias{downloadIfOld}
\title{Download a file}
\usage{
downloadIfOld(url, file = basename(url), age = "5 days", verbose = FALSE,
  source = FALSE)
}
\arguments{
\item{url}{a string specifying the remote location}

\item{file}{local file name}

\item{age}{maximum age of the local file}

\item{verbose}{print additional information}

\item{source}{source the file after download}
}
\description{
Checks if a downloaded file is old and re-downloads is necessary
}

