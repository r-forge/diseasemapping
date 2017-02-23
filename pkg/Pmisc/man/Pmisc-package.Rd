\name{Pmisc-package}
\alias{Pmisc-package}
\alias{Pmisc}
\docType{package}
\title{
Pmisc package
}
\description{
Various utilities which will eventually be put into other packages
}


\author{
Patrick Brown
}

\examples{


mydir = tempdir()

file.copy(
	system.file("docs/document.Rmd", package='Pmisc'),
	mydir, overwrite = TRUE)

Makefile("document", suffix=c('tmd', 'pdf','html','docx'),
	output = file.path(mydir, 'Makefile'), 
	cygwin = .Platform$OS.type == 'windows')

# if necessary, reformat directory for cygwin
mydirMake = gsub("\\\\", "/", gsub("^C:", "/cygdrive/c", mydir))

\dontrun{

# run make
system(
  paste("bash -c \"cd ", mydirMake, "; make\"", sep="") 
)
	
}

file.path(mydir, 'document.pdf')

}
