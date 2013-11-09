\name{colourScale}
\alias{colourScale}

\title{Create colour scales}
\description{
Produces a scale of colours for plotting maps an producing legends.  }
\usage{
	colourScale(x, breaks=5, style=c("quantile","equal","unique", "fixed"),
	 col="YlOrRd", opacity=1, dec=NULL, firstBreak=NULL, 
	 transform=NULL, revCol=FALSE, exclude=NULL)
	}

\arguments{
	\item{x}{A vector or single-layer Raster, numeric or factor,  for which a colour scale will be created}
	\item{breaks}{Either the number of breaks or a vector of breaks}
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
	\item{...}{Additional arguments passed to \code{\link[classInt]{classIntervals}}}
}

\details{
Produces intervals from \code{x}, each with a unique colour. Categories are determined with break points according to the following \code{style} options:
\itemize{
\item \code{quantile}: \code{quantile(x, prob=seq(0,1,len=breaks}, )} 
\item \code{equal}: \code{seq(min(x), max(x), len=breaks)} 
\item \code{unique}: \code{sort(table(unique(x)))[1:breaks]} 
\item \code{fixed}: \code{breaks}
\item any other string: is passed to \code{\link[classInt]{classIntervals}}
}
}

\value{A list with elements
\item{col}{Vector of same length of \code{x} containing colours (RGB hex)}
\item{breaks}{vector of break points}
\item{legendCol}{vector of unique colour values corrdsponding to \code{breaks}}
}
\seealso{\code{\link{legendBreaks}}, \code{\link[classInt]{classIntervals}}

\examples{
Npoints = 20
myPoints = SpatialPointsDataFrame(20*cbind(runif(Npoints), runif(Npoints)), 
	data=data.frame(y1=c(NA, rnorm(Npoints-1)), 
	y2=c(sample(1:4, Npoints-1,replace=TRUE), NA)),
	proj4string=CRS("+init=epsg:4326")
	)
mymap = openmap(myPoints)


myscale = colourScale(myPoints$y1, breaks=4, col="RdYlBu", style="quantile", revCol=TRUE,dec=1)



map.new(myPoints)
plot(mymap,add=TRUE)
plot(myPoints, col=myscale$col, pch=16,add=TRUE)
legendBreaks("topleft", breaks=myscale$breaks, col=myscale$legendCol)


myscale2 = colourScale(myPoints$y1, breaks=8, col=rainbow, style="equal",
opacity=0.8, dec=2, revCol=TRUE)

map.new(myPoints)
plot(mymap,add=TRUE)
plot(myPoints, col=myscale2$col, pch=16,add=TRUE)
legendBreaks("topleft", breaks=myscale2$breaks, col=myscale2$legendCol)


myscale3 = colourScale(myPoints$y2,  breaks=3,col="Set2", style="unique",
	opacity=c(0.1, 0.9))

map.new(myPoints)
plot(mymap,add=TRUE)
plot(myPoints, col=myscale3$col, pch=16,add=TRUE)
legendBreaks("topleft", breaks=myscale3$breaks, col=myscale3$legendCol)

myPoints$y3 = exp(myPoints$y1)
myscale4 = colourScale(myPoints$y3,  breaks=4, style="equal",
	opacity=c(0.1, 0.9), transform=1.25,dec=3)

map.new(myPoints)
plot(mymap,add=TRUE)
plot(myPoints, col=myscale4$col, pch=16,add=TRUE)
legendBreaks("topleft", breaks=myscale4$breaks, col=myscale4$legendCol)

}