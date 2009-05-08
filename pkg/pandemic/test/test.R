source("../R/parameters.R")
source("../R/simEpidemic.R")

params = pandemicParams()
 data = simEpidemic(params, 10)
 summary(data)

 xseq=0:20

plot(xseq, dweibullRound(xseq, params$HospRec))
sum( dweibullRound(0:100, params$HospRec))

priors = pandemicPriors()
names(priors)
names(priors$InfOns)
priors$InfOns$mean
