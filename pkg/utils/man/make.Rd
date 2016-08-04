\name{make}
\alias{make}
\title{Run make}
\usage{
make(x = commandArgs(TRUE), suffix = NULL, citeproc = 0,
  beamer = FALSE, run = FALSE, ...)
}
\arguments{
\item{x}{a string specifying the target}

\item{suffix}{string to replace file extension}

\item{citeproc}{use \code{--filter=pandoc-citeproc}}

\item{beamer}{use \code{--to=beamer}}

\item{run}{execute the command}

\item{...}{additional arguments passed to make}
}
\description{
Runs make with the Makefile supplied by \code{Pmisc}
}

\details{

The following can be used in an eclipse configuration

\preformatted{--vanilla -e "Pmisc::make(suffix='html', run=TRUE)" --args ${selected_resource_name}}

}

\examples{
make("stuff.html")
}
