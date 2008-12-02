\name{mergeBugsData}
\alias{mergeBugsData}
\alias{mergeBugsData.SpatialPolygonsDataFrame}
\alias{mergeBugsData.data.frame}


\title{ merge the result from bugs function}
\description{
 merge the result from bugs function
}
\usage{
mergeBugsData(x, bugsSummary, by.x = NULL, newcol = "mean", ...)
\method{mergeBugsData}{SpatialPolygonsDataFrame}(x, bugsSummary, by.x=NULL, newcol="mean", ...)
\method{mergeBugsData}{data.frame}(x, bugsSummary, by.x=NULL, newcol="mean", ...)
}

\arguments{
  \item{x}{spatial polygon object i.e population data set (popdata)}
  \item{bugsSummary}{ posterior distribution result from summaryChain function }
  \item{by.x}{ the common term from the spatial polygon object and the bugs function result  }
  \item{newcol}{the summay statistic that to be merged back to the data frame  }
  \item{\dots}{ additional arguments }
}

\author{ Patrick Brown }


\examples{
data(popdata)
newdata = c("3560102"=2, "3560104"=3)
popdatatry = mergeBugsData(popdata, newdata, by.x="CSDUID")

# if the data set is a spatial polygons data frame: 
#popdatatry = mergeBugsData.SpatialPolygonsDataFrame(popdata, newdata, by.x="CSDUID")

# if the data set is a data frame
#popdatatry = mergeBugsData.data.frame(popdata, newdata, by.x="CSDUID")

\dontrun{
library(glmmBUGS)
data(popdata)
data(casedata)

model = getRates(casedata, popdata, ~age*sex)
ontario = getSMR(model, popdata, casedata)
ontario = ontario@data[,c("CSDUID","Observed","logExpected")]
popDataAdjMat = poly2nb(ontario,row.names=as.character(ontario[["CDSUID"]]))
data(popDataAdjMat)

forBugs = glmmBUGS(formula=Observed + logExpected ~ 1,
  effects="CSDUID",   family="poisson", spatial=popDataAdjMat,
  data=ontario)

startingValues = forBugs$startingValues

source("getInits.R")

ontarioResult = bugs(forBugs$ragged, getInits, parameters.to.save = names(getInits()),
    model.file="model.bug", n.chain=3, n.iter=100, n.burnin=10, n.thin=2,
      program="winbugs", debug=TRUE)

data(ontarioResult)
ontarioParams = restoreParams(ontarioResult, forBugs$ragged)
ontarioSummary = summaryChain(ontarioParams)
ontario = mergeBugsData(ontario, ontarioSummary)
}

}
