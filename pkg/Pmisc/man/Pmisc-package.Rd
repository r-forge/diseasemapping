\name{Pmisc-package}
\alias{Pmisc-package}
\alias{Pmisc}
\docType{package}
\title{
Pmisc package
}
\description{
Various functions for knitr documents, INLA, and lgcp.
}

\details{

The example below makes html, docs, and pdf files from a simple template \code{document.Rmd}.

}

\author{
Patrick Brown
}

\examples{
library("Pmisc")

list.files(system.file("docs", package="Pmisc"), full.names=TRUE)
file.show(system.file("docs/document.Rmd", package="Pmisc"))
# copy the document.Rmd file to the temporary folder
(mydir = tempdir())
file.copy(
	system.file("docs/document.Rmd", package='Pmisc'),
	mydir, overwrite = TRUE)

# create a Makefile for building document.Rmd into pdf, html and docx
Makefile("document", suffix=c('pdf','html','docx'),
	output = file.path(mydir, 'Makefile'))


# run make to build the output files
if(interactive()) {
# if necessary, reformat directory name for cygwin
	mydirMake = gsub("^C:", "/cygdrive/c", normalizePath(mydir, "/"))
	commandString = paste("cd ", mydirMake, "; make")
	system2("bash", paste('-c "', commandString, '"'))
}

# show the files
list.files(mydir, 
  '*.(pdf|docx|html|.md|tex)$',
  full.names=TRUE, ignore.case=TRUE)  
}
