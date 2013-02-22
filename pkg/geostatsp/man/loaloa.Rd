\name{loaloa}
\alias{loaloa}
\alias{elevatinLoa}
\alias{eviLoa}
\docType{data}
\title{
Loaloa prevalence data from 197 village surveys
}
\description{
Location and prevalence data from villages, elevation an vegetation index for the study region.
}
\usage{data(loaloa)}
\format{
\code{loaloa} is a SpatialPolygonsDataFrame of the data, with columns \code{N} being the number 
of individuals tested and \code{y} being the number of positives.
\code{elevationLoa} is a raster of elevation data.
\code{elevationLoa} is a raster of vegetation index.
}
\source{
\url{http://www.leg.ufpr.br/doku.php/pessoais:paulojus:mbgbook:datasets} for the loaloa data,
\url{http://e4ftl01.cr.usgs.gov/MOLT/MOD13Q1.005/} for the EVI data,  and  
\url{http://srtm.csi.cgiar.org} for the elevaiton data.
}

\examples{
data(loaloa)
plot(eviLoa)
points(loaloa)





}


\keyword{datasets}
