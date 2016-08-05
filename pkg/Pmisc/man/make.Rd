\name{make}
\alias{make}
\title{Run make}
\usage{
make(x, suffix = NULL,
  beamer = FALSE, run = FALSE, ...)
}
\arguments{
\item{x}{a string specifying the target}

\item{suffix}{string to replace file extension}

\item{beamer}{use \code{--to=beamer}}

\item{run}{execute the command}

\item{...}{additional arguments passed to make}
}
\description{
Runs make with the Makefile supplied by \code{Pmisc}
}

\details{

The following can be used in an eclipse configuration

\preformatted{\bin\bash -c "make $(R --vanilla --slave -e Pmisc::make\(\'${selected_resource_name}\',suffix=\'pdf\'\))"}

}

\examples{
make("stuff.html")

make("stuff.pdf", beamer=TRUE)

}
