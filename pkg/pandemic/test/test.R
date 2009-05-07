source("../R/parameters.R")
source("../R/simEpidemic.R")

params = pandemicParams()
 data = simEpidemic(params, 10)
 summary(data)
