\name{legendBreaks}
\alias{legendBreaks}

\title{Legends for colour scale}
\description{
Legends where N+1 labels are supplied as the limits of N bins.
}
\usage{
legendBreaks(pos, 
	breaks,     
	col=breaks$col,
    legend=breaks$breaks,
    rev=TRUE,
    outer=TRUE,
    pch=15,
    cex=par('cex'),
    pt.cex=2*cex,
    adj=c(0,0.75)*par('cra')*cex,
    y.intersp=0.75,
    inset=0.05,
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
\item{cex}{see \code{\link[graphics]{legend}}}
\item{pt.cex}{see \code{\link[graphics]{legend}}}
\item{adj}{see \code{\link[graphics]{legend}}}
\item{y.intersp}{see \code{\link[graphics]{legend}}}
\item{inset}{see \code{\link[graphics]{legend}}}
\item{...}{Additional arguments passed to \code{\link[graphics]{legend}}.}
}

\details{


A legend for 'z-axis' colour scales.


}

\value{
Result of call to \code{\link[graphics]{legend}}
}

\seealso{\code{\link{colourScale}}}

