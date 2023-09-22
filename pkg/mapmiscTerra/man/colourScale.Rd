\name{colourScale}
\alias{colourScale}
\alias{colorScale}
\alias{breaksForRates}

\title{Create colour scales}
\description{
Produces a scale of colours for plotting maps}
\usage{
	colourScale(x, breaks=5, style=c("quantile","equal","unique", "fixed"),
	 col="YlOrRd", opacity=1, dec=NULL, digits = 6, firstBreak=NULL, 
	 transform=NULL, revCol=FALSE, exclude=NULL, labels=NULL, ...)
	colorScale(...)
	breaksForRates(x, breaks = 10, transform = 0.1,
      multiples = c(2, 4, 5, 10)) 
	}

\arguments{
	\item{x}{A vector or single-layer Raster, numeric or factor,  for which a colour scale will be created}
	\item{breaks}{For \code{colourScale} either the number of or vector of breaks.
	for \code{legendBreaks} usually the output of \code{colourScale}, or a vector of breaks}
	\item{style}{Style for breaks, see Details}
	\item{col}{Colours to use, either a function or 
	argument for \code{\link[RColorBrewer:ColorBrewer]{brewer.pal}}}
	\item{opacity}{adds transparency to colours, either a single number,
	vector of length 2, or vector of same length as \code{breaks}}
	\item{dec}{Number of decimal places for the breaks}	
	\item{digits}{Number of significant figures}
	\item{firstBreak}{If non-null, force the first break to take this value (often zero).}
	\item{transform}{A list of two functions to transform \code{x} and inverse transform the breaks, or a numeric value specifying a Box-Cox parameter.}
	\item{revCol}{Reverse the order of the colours.}
	\item{exclude}{A vector of values to change to NA when they appear in \code{x}}
	\item{labels}{Vector of names of levels, useful when \code{style='unique'}}
	\item{multiples}{break points must be multiples of these numbers times a power of 10}
	\item{...}{Additional arguments passed to \code{\link[classInt]{classIntervals}}.}
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

\code{colorScale} passes all it's arguments to \code{colourScale}

\code{breaksForRates} returns break points suitable for mapping incidence rates, which are positive and always
include 1.0.

}

\value{A list with elements
\item{plot}{Vector of same length of \code{x} containing colours (RGB hex)}
\item{breaks}{vector of break points}
\item{col}{vector of unique colour values corresponding to \code{breaks}}
\item{colWithOpacity}{as \code{col}, but with two digit transparency values appended.}
}

\seealso{\code{\link{legendBreaks}},\code{\link{scaleBar}}, \code{\link[classInt]{classIntervals}}}

\examples{

breaksForRates(13.6, breaks = 7) 

Npoints = 20

myPoints = vect(
	cbind(runif(Npoints), 51+runif(Npoints)), 
	atts=data.frame(
		y1=c(NA, rnorm(Npoints-1)), 
	    y2=c(sample(0:5, Npoints-1,replace=TRUE), NA)
	),
	crs=crsLL)



if(require('RColorBrewer', quietly=TRUE)) {
	theCol = 'RdYlBu'
} else {
	theCol = grDevices::heat.colors
}

myscale = colourScale(myPoints$y1, breaks=4, col=theCol,
 style="quantile", revCol=TRUE,dec=1)


data("netherlands")
nldElev = terra::unwrap(nldElev)
myscale = colourScale(nldElev, breaks=4, col=theCol, style='equal', dec=0)

oldpar = map.new(myPoints)
plot(myPoints, col=myscale$plot, pch=16,add=TRUE)
legendBreaks("topleft", myscale)


myscale2 = colourScale(myPoints$y1, breaks=8, col=rainbow, style="equal",
opacity=0.8, dec=2, revCol=TRUE)

map.new(myPoints)
plot(myPoints, col=myscale2$plot, pch=16,add=TRUE)
legendBreaks("topleft", myscale2)

if(require('RColorBrewer', quietly=TRUE)) {
	theCol = 'Set2'
} else {
	theCol = grDevices::heat.colors
}

myscale3 = colourScale(myPoints$y2,  breaks=3,col=theCol, style="unique",
	opacity=c(0.1, 0.9))

map.new(myPoints)
plot(myPoints, col=myscale3$plot, pch=16,add=TRUE)
legendBreaks("topleft", myscale3)

myPoints$y3 = exp(myPoints$y1)
myscale4 = colourScale(myPoints$y3,  breaks=4, style="equal",
	opacity=c(0.1, 0.9), transform=1.25,dec=0, firstBreak=0)

map.new(myPoints)
plot(myPoints, col=myscale4$plot, pch=16,add=TRUE)
legendBreaks("topleft", myscale4)


# raster with colour table

x = rast(extent=ext(0,15,0,10), res=1)
values(x) = sample(1:4, ncell(x), replace=TRUE)
myScale = colourScale(x, breaks=3, style='unique', col=c('red','blue','orange'))
if(utils::packageVersion("terra") >= "1.7-40" ) {
	terra::coltab(x) = myScale$colourtable
	plot(x)
} else {
	plot(x, breaks = myScale$breaks, col=myScale$col)
}
legendBreaks('topright', myScale)

par(oldpar)
}





