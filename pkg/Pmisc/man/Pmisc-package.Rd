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


\author{
Patrick Brown
}

\examples{
library("Pmisc")

# copy the document.Rmd file to the temporary folder
(mydir = tempdir())
file.copy(
	system.file("docs/document.Rmd", package='Pmisc'),
	mydir, overwrite = TRUE)

# create a Makefile for building document.Rmd into pdf, html and docx
Makefile("document", suffix=c('tmd', 'pdf','html','docx'),
	output = file.path(mydir, 'Makefile'), 
	cygwin = .Platform$OS.type == 'windows')


# run make to build the output files
if(interactive()) {
# if necessary, reformat directory name for cygwin
	mydirMake = gsub("\\\\", "/", 
	  gsub("^C:", "/cygdrive/c", mydir))
	system(paste(
	  "bash -c \"cd ", mydirMake, 
	  "; make\"", sep=""))
}

# show the files
list.files(mydir, 
  '*.(pdf|docx|html|.md|tex)$',
  full.names=TRUE, ignore.case=TRUE)  
}
