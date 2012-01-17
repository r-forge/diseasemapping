\name{bym}
\alias{bym}

\title{Fit a Besag-York-Mollie Model using INLA}
\arguments{
  \item{formula}{the formula used to specify the model, default is "observed ~ 1", intercept only model}
  \item{family}{Distribution of responses, default is "poisson" }
  \item{spPolygon}{name of spatial polygon, with data attribute, containing observed and expected case per region}
  \item{offset}{offest variable name in the dataset, default is "logExpected"}
  \item{ID}{Identification column of the regions, which follows GMRF, default is "ID" }
  \item{prior.Spatial}{A string specify the prior distribution  for spatail random effect,
                          default is "c(1,0.01)"}
  \item{prior.iid}{A string specify the prior distribution for iid effect,
                          default is "c(1,0.01)"}
  \item{exceed}{A number bigger than 1 for exceedance probability, default is 1.5, compute 50% exceedance probability}
  \item{debug}{Print Debug information, default is FALSE}
  \item{verbose}{Print program excution process, default is TRUE}
  \item{keep}{Keep the data and result files, default is TRUE}
}
\details{
It calls INLA packages, obtain the marginal distribution of the spatial random effect and the iid+spatial random effects
}
\value{
 A list with a spatial polygon object and an INLA object
}
\references{  }
\author{ Patrick Brown  }
\note{
}
\seealso{ }
\examples{
 data(popdata)
 data(casedata)
 model = getRates(casedata, popdata, ~age*sex, breaks=seq(0, 90, by=10) )
 ontario = getSMR(popdata,model,  casedata,regionCode="CSDUID", regionCodeCases="CSD2006")
result =bym(formula=observed  ~ 1,ID="CSDUID", family="poisson", spPolygon=ontario@data,offset = "logExpected")
)



###Run BYM on grid
library(spdep)
nb8<-cell2nb(nrow=50, ncol=50, type="queen", torus=T)
nb2inla("germany.graph",nb8) #Write graph file, then run bym function
}

\keyword{INLA}
\keyword{bym}% __ONLY ONE__ keyword per line
