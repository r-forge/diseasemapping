\name{colourScale}
\alias{colourScale}
\alias{legendBreaks}

\title{Create colour scales}
\description{
Produces a scale of colours for plotting maps an producing legends.  }
\usage{
	colourScale(x, breaks=5, style=c("quantile","equal","unique", "fixed"),
	 col="YlOrRd", opacity=1, dec=NULL, firstBreak=NULL, 
	 transform=NULL, revCol=FALSE, exclude=NULL, ...)
legendBreaks(pos, breaks, outer=TRUE,...)
	}

\arguments{
	\item{x}{A vector or single-layer Raster, numeric or factor,  for which a colour scale will be created}
	\item{breaks}{For \code{colourScale} either the number of or vector of breaks.
	for \code{legendBreaks} usually the output of \code{colourScale}, or a vector of breaks}
	\item{style}{Style for breaks, see Details}
	\item{col}{Colours to use, either a function or 
	argument for \code{\link[RColorBrewer]{brewer.pal}}}
	\item{opacity}{adds transparency to colours, either a single number,
	vector of length 2, or vector of same length as \code{breaks}}
	\item{dec}{Number of decimal places for the breaks}
	\item{firstBreak}{If non-null, force the first break to take this value (often zero).}
	\item{transform}{A list of two functions to transform \code{x} and inverse transform the breaks, or a numeric value specifying a Box-Cox parameter.}
	\item{revCol}{Reverse the order of the colours.}
	\item{exclude}{A vector of values to change to NA when they appear in \code{x}}
  \item{pos}{
Position, as specified in the \code{legend} function.
}
\item{outer}{If \code{TRUE}, put legend in the margin of the plot}
	\item{...}{Additional arguments passed to \code{\link[classInt]{classIntervals}} or \code{legend}.}
}

\details{
\code{colourScale} produces intervals from \code{x}, each with a unique colour. Categories are determined with break points according to the following \code{style} options:
\itemize{
\item \code{quantile}: \code{quantile(x, prob=seq(0,1,len=breaks), )} 
}

 \code{equal}: \code{seq(min(x), max(x), len=breaks)} 


 \code{unique}: \code{sort(table(unique(x)))[1:breaks]} 


 \code{fixed}: \code{breaks}


 any other string: is passed to \code{\link[classInt]{classIntervals}}

\code{legendBreaks} produces a legend for 'z-axis' colour scales.


}

\value{A list with elements
\item{plot}{Vector of same length of \code{x} containing colours (RGB hex)}
\item{breaks}{vector of break points}
\item{col}{vector of unique colour values corresponding to \code{breaks}}
\item{colWithOpacity}{as \code{col}, but with two digit transparency values appended.}
}

\seealso{\code{\link{scaleBar}}, \code{\link[classInt]{classIntervals}}}

\examples{
Npoints = 20
myPoints = SpatialPointsDataFrame(20*cbind(runif(Npoints), runif(Npoints)), 
	data=data.frame(y1=c(NA, rnorm(Npoints-1)), 
	y2=c(sample(1:4, Npoints-1,replace=TRUE), NA)),
	proj4string=CRS("+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0")
	)

\dontrun{
mymap = openmap(myPoints)
}

if(require('RColorBrewer', quietly=TRUE)) {
	theCol = 'RdYlBu'
} else {
	theCol = heat.colors
}

myscale = colourScale(myPoints$y1, breaks=4, col=theCol,
 style="quantile", revCol=TRUE,dec=1)



map.new(myPoints)
\dontrun{
plot(mymap,add=TRUE)
}
plot(myPoints, col=myscale$plot, pch=16,add=TRUE)
legendBreaks("topleft", breaks=myscale)


myscale2 = colourScale(myPoints$y1, breaks=8, col=rainbow, style="equal",
opacity=0.8, dec=2, revCol=TRUE)

map.new(myPoints)
\dontrun{
plot(mymap,add=TRUE)
}
plot(myPoints, col=myscale2$plot, pch=16,add=TRUE)
legendBreaks("topleft", breaks=myscale2)

if(require('RColorBrewer', quietly=TRUE)) {
	theCol = 'Set2'
} else {
	theCol = heat.colors
}

myscale3 = colourScale(myPoints$y2,  breaks=3,col=theCol, style="unique",
	opacity=c(0.1, 0.9))

map.new(myPoints)
\dontrun{
plot(mymap,add=TRUE)
}
plot(myPoints, col=myscale3$plot, pch=16,add=TRUE)
legendBreaks("topleft", breaks=myscale3)

myPoints$y3 = exp(myPoints$y1)
myscale4 = colourScale(myPoints$y3,  breaks=4, style="equal",
	opacity=c(0.1, 0.9), transform=1.25,dec=0, firstBreak=0)

map.new(myPoints)
\dontrun{
plot(mymap,add=TRUE)
}
plot(myPoints, col=myscale4$plot, pch=16,add=TRUE)
legendBreaks("topleft", breaks=myscale4$breaks, col=myscale4$col)

}