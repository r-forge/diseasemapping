\name{markdownHeader}
\alias{markdownHeader}
\alias{today}
\title{Markdown headers}
\usage{
markdownHeader(title = NULL, author = NULL, date = Pmisc::today(),
  startAndEnd = "---", slides = TRUE, ...)
today(...)
}
\arguments{
\item{title}{value for title}

\item{author}{value for author field}

\item{date}{value for date field}

\item{startAndEnd}{Add \code{---} before and after YAML block}

\item{slides}{Add commands intended for producing slides.}

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
cat(markdownHeader(
	title='my title',
	author = 'me',
	extrastuff = c(a=1, b=2)
))
}

