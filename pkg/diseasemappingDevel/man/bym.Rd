\name{bym}
\alias{bym}

\title{Fit a Besag-York-Mollie Model using INLA}
\arguments{
  \item{formula}{the formula used to specify the model, default is "observed ~ 1", intercept only model}
  \item{family}{Distribution of responses, default is "poisson" }
  \item{data}{dataset name, containing observed and expected case per region}
  \item{offset}{offest variable name in the dataset, default is "logExpected"}
  \item{spatial}{Identification column of the regions, which follows GMRF, default is "region" }
  \item{graph.file}{A file contains the neigbour structures, see INAL package}
  \item{control.Spatial}{A string specify the prior distribution and starting values for spatail random effect,
                          default is "param=c(1,0.01),initial=4"}
  \item{control.iid}{A string specify the prior distribution and starting values for iid effect,
                          default is "param=c(1,0.01),initial=4"}
  \item{debug}{Print Debug information, default is FALSE}
  \item{verbose}{Print program excution process, default is TRUE}
  \item{keep}{Keep the data and result files, default is TRUE}
}
\details{
It calls INLA packages, obtain the joint distrbution of spatail term and non-spatail term
using linear combinations.
}
\value{
 A INLA object
}
\references{  }
\author{ Patrick Brown  }
\note{
}
\seealso{ }
\examples{
library(INLA)
data(Germany)
Germany$logE = log(Germany$E)
a<-bym(formula="Y ~ x",spatial="region",offset = "logE",graph.file="germany.graph",data=Germany)
a$summary.fixed


#To merge the result back to the data. This will add  7 additional columns in the dataset
#Posterior means of spatial effects, iid effects, fitted values, fitted values on the log scale,
#bym (the spatial and non-spatial effects), as well as the exceedance probability of fitted and bym values
resultData<-mergeINLA(data=Germany,INLA=a,exceed=1.2)

###Run BYM on grid
library(spdep)
nb8<-cell2nb(nrow=50, ncol=50, type="queen", torus=T)
nb2inla("germany.graph",nb8) #Write graph file, then run bym function
}

\keyword{INLA}
\keyword{bym}% __ONLY ONE__ keyword per line
