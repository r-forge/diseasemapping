\name{markdownHeader}
\alias{markdownHeader}
\alias{today}
\alias{mathCommands}
\alias{subcaptionCommands}
\title{Markdown headers}
\usage{
markdownHeader(
	title=NULL, 
	author=NULL,
	date=Pmisc::today(),
	biblatex=1,
	bibliotitle = 'References',
	bibliostyle = 'authoryear,backend=biber',
	biblatexoptions = c(
		maxbibnames=20,
		maxcitenames=2,doi='false',
		isbn='false',	
		giveninits='true',
		uniquelist='false'),
	beamer=FALSE,
	mathCommands=FALSE,
	subcaptionCommands = !beamer,
	...
)
today(...)
mathCommands
subcaptionCommands
}
\arguments{
\item{title}{value for title}
\item{author}{value for author field}
\item{date}{value for date field}
\item{biblatex}{uses the biblatex package in LaTeX}
\item{bibliotitle}{title for references section}
\item{bibliostyle}{style of bibliography}
\item{biblatexoptions}{options for citations, passed to the biblatex package}
\item{beamer}{Add commands intended for producing slides with LaTeX beamer.}
\item{mathCommands}{Defines latex commands for var, E, tp, dens, and a few others.}
\item{subcaptionCommands}{Enables use of the subcaption package with knitr.}
\item{...}{Additional items for the header, such as 
\code{header-includes}, or for \code{\link[base]{format}} 
when passed to \code{today}.
}
}
\value{
A string containing the header.  \code{today()} returns today's date
}
\description{
Creates YAML headers for Markdown slides via pandoc
}
\details{
Put this function in a code chunk at the start of
}
\examples{
today()
cat(Pmisc::markdownHeader(
	title='my title',
	author = 'me',
	geometry='margin=1in',
	fontsize='12pt',
	headerincludes=c('\\newcommand{\\tp}{^{\\textrm{T}}}')
))
}

