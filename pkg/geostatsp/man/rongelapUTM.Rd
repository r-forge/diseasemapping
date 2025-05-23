\name{rongelapUTM}
\alias{rongelapUTM}
\docType{data}
\title{
Rongelap data
}
\description{
  This data-set was used by Diggle, Tawn and Moyeed (1998) to illustrate
  the model-based geostatistical methodology introduced in the paper.
  discussed in the paper. The radionuclide concentration data set consists 
  of measurements of \eqn{\gamma}{gamma}-ray counts at \eqn{157} locations. 
}
\usage{data(rongelapUTM)}
\format{
A \code{SpatVector}, with columns \code{count} being the 
radiation count and \code{time} being the length of time
the measurement was taken for.  A UTM coordinate reference system is used, where coordinates are in metres.
}

\source{
\url{https://web.archive.org/web/20240110054727/http://www.leg.ufpr.br/doku.php/pessoais:paulojus:mbgbook:datasets}.  For further details on the radionuclide concentration data, see Diggle,Harper and Simon (1997), Diggle, Tawn and Moyeed (1998) and Christensen (2004).
}

\examples{
data("rongelapUTM")
rongelapUTM = unwrap(rongelapUTM)
plot(rongelapUTM, main="Rongelap island")

if(require('mapmisc')) {
  bgMap = openmap(rongelapUTM, buffer=300, maxTiles=2)
  plot(bgMap)
  points(rongelapUTM, cex=0.4)
  scaleBar(rongelapUTM, 'left')

}



}

\references{
  Christensen, O. F. (2004). Monte Carlo maximum likelihood in model-based geostatistics.
  Journal of computational and graphical statistics \bold{13} 702-718.
  
  Diggle, P. J., Harper, L. and Simon, S. L. (1997). Geostatistical analysis of residual contamination from nuclea
  testing. In: \emph{Statistics for the environment 3: pollution assesment and control} (eds. V. Barnet and K. F. Turkmann), Wiley, 
  Chichester, 89-107.
  
  Diggle, P. J., Tawn, J. A. and Moyeed, R. A. (1998). Model-based
  geostatistics (with Discussion). Applied Statistics,
  47, 299--350.
}


\keyword{datasets}
