\name{col2html}
\alias{col2html}

\title{Convert colours to HTML hex}
\description{
Converts any object interpretable as a colour to an HTML hex string, i.e. 'red' to '#FF0000'.
}
\usage{
col2html(col, opacity=1, alpha)
	}

\arguments{
\item{col}{Either a character vector of colour names 
as listed by \code{\link[grDevices:colors]{colours}()}
 or an integer vector of colour indexes. Passed to \code{\link[grDevices]{col2rgb}}.}
\item{opacity}{scalar or vector of colour opacities between 0 and 1.}
\item{alpha}{Integer between 0 and 255, or a character giving a 2-digit hex value.  
Overrides \code{opacity} and passed to \code{\link[grDevices]{rgb}}.}
}

\value{
A vector of 6 or 8 digit hex codes specifying HTML colours.
}

\seealso{\code{\link[grDevices]{col2rgb}},\code{\link[grDevices]{rgb}}\code{\link[base]{hexmode}}}

\examples{
col2html(1:10)
col2html(c('red','blue'),0.5)
col2html(c(2,4),0.5)
col2html(c(stuff='red',foo='blue'),alpha=128)
col2html(c('red','blue'),alpha='80')
col2html(c(2,4),alpha='80')

N = length(palette())
plot(1:N, rep(1,N),xlim=c(0,N),pch=16,cex=5,
col=col2html(1:N))
points(1:N, rep(1,N),pch=15,cex=4.5, col=palette())
text(-0.5+1:10, rep(1,10), col2html(1:10),srt=90)
text(1:N, rep(0.7,N), palette())
text(1:N-0.5, rep(1.3, N), col2html(palette()), cex=0.7)
}