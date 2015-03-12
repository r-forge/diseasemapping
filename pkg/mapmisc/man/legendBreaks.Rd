\name{legendBreaks}
\alias{legendBreaks}

\title{Legends for colour scale}
\description{
Legends where N+1 labels are supplied as the limits of N bins.
}
\usage{
legendBreaks(pos, 
	breaks,     
	col,    legend,
    rev=TRUE,
    outer=TRUE,
    pch=15,
    bg=par('bg'),
    cex=par('cex'),
    pt.cex=2.5*cex,
    text.col=par('fg'),
    title=NULL,
    inset=0.05,
    title.col=text.col,
    adj=0,
    y.intersp,
	...)
	}

\arguments{
  \item{pos}{
Position, as specified in the \code{\link[graphics]{legend}} function.
}
	\item{breaks}{Optional list with elements \code{col} and \code{legend}, such as the 
	output from  \code{\link{colourScale}} }
	\item{col}{Single colour or vector of colours for each bin}
	\item{legend}{vector of labels for the legend, one more element than there are 
	colours}
\item{rev}{if \code{TRUE}, labels and colours are ordered from bottom to top, otherwise top to bottom. }
\item{outer}{If \code{TRUE}, put legend in the margin of the plot}
\item{pch}{see \code{\link[graphics]{legend}}}
\item{bg}{background colour see \code{\link[graphics]{legend}}}
\item{cex}{see \code{\link[graphics]{legend}}}
\item{pt.cex}{see \code{\link[graphics]{legend}}}
\item{text.col}{see \code{\link[graphics]{legend}}}
\item{title}{see \code{\link[graphics]{legend}}}
\item{inset}{see \code{\link[graphics]{legend}}}
\item{title.col}{see \code{\link[graphics]{legend}}}
\item{adj}{Adjustment of the legend labels relative to plotting symbols.}
\item{y.intersp}{see \code{\link[graphics]{legend}}}
\item{...}{Additional arguments passed to \code{\link[graphics]{legend}}.}
}

\details{


A legend for 'z-axis' colour scales.


}

\value{
Result of call to \code{\link[graphics]{legend}}
}

\seealso{\code{\link{colourScale}}}

