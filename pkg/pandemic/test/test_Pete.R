source("../R/parameters.R")
source("../R/simEpidemic.R")
source("../R/weibullRound.R")

params = pandemicParams()
 xdata = simEpidemic(params, 100)
data=xdata

source("../R/dataInitialize.R")

source("../R/mcmcCode.R")

data

 