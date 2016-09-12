\name{make}
\alias{make}
\alias{Makefile}
\title{Create a Makefile}
\usage{
make(x, suffix = NULL,
  beamer = FALSE, run = FALSE, ...)
Makefile(x, suffix = NULL,
  beamer = FALSE, output='Makefile', cygwin=FALSE)
}
\arguments{
\item{x}{a string specifying the target}

\item{suffix}{string to replace file extension}

\item{beamer}{use \code{--to=beamer}}

\item{run}{run make with a system call}

\item{output}{Name of file saved, or a function output is passed to.}

\item{cygwin}{Makefile is to be run in cygwin}

\item{...}{additional arguments passed to make}
}
\description{
Creates a Makefile or a string to run \code{make} with the Makefile supplied by \code{Pmisc}
}

\details{

The following can be used in an eclipse configuration

\preformatted{\bin\bash -c make $(}
\preformatted{"R --vanilla --slave -e }
\preformatted{Pmisc::make\(\'${selected_resource_name}\',suffix=\'pdf\'}
\preformatted{\))"}

}

\examples{
Makefile(x='stuff.md', suffix='pdf', beamer=TRUE, output=cat)

make("stuff.html")

stuff = make("stuff.pdf", beamer=TRUE)

stuff

if(.Platform$OS.type == 'unix') {
	system(paste(
	"echo make $(R --vanilla --slave -e",
	"\"Pmisc::make('stuff.Rmd',suffix='pdf')\")")
	)
}

}
