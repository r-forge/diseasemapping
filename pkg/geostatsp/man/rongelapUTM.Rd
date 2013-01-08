\name{rongelapUTM}
\alias{rongelapUTM}
\docType{data}
\title{
Rongelap data
}
\description{
Data on radiation measurements on Rongelap island, as a SpatialPointsDataFrame in 
a UTM coordinate reference system.
}
\usage{data(rongelapUTM)}
\format{
A SpatialPolygonsDataFrame, with columns \code{count} being the 
radiation count and \code{time} being the length of time
the measurement was taken for.
}

\source{
http://www.leg.ufpr.br/doku.php/pessoais:paulojus:mbgbook:datasets
}

\examples{
data(rongelapUTM)
plot(rongelapUTM)
\dontrun{
load(url("http://www.filefactory.com/file/54e8egxfddul/n/MHL_adm0_RData"))

rongelapLL = spTransform(rongelapUTM, gadm@proj4string)
plot(rongelapLL)
plot(gadm, add=T)
}


}
\keyword{datasets}
