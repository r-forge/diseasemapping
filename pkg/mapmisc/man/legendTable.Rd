\name{legendTable}
\alias{legendTable}

\title{Table for colour scales}
\description{
A table in html or Latex showing values associated with colours
}
\usage{
legendTable(x,
    type=c('latex',  'html'),
    box = c(-0.2, 1, 2),
    unit = 'em',
    collapse=NULL)
}

\arguments{
  \item{x}{
a \code{data.frame} with columns \code{col} and \code{label}, possibly 
produced by \code{\link{colourScale}}
}
\item{type}{html or latex compatible output}
\item{box}{dimensions of colour boxes, passed as 
depth, height and width to \code{rule} in Latex, or 
width (first two elements ignored) for html. 
}
\item{unit}{Units for box dimensions}
\item{collapse}{If non-NULL, passed to \code{paste} to produce a character 
vector instead of table}
}



\value{
\code{data.frame} or character vector
}

\seealso{\code{\link{colourScale}}}

\examples{

mytable = data.frame(col=col2html(1:5), label=1:5)

legendTable(mytable)
legendTable(mytable, collapse=';')
legendTable(mytable, type='html')

}

