source("../R/parameters.R")
source("../R/simEpidemicInfection.R")
source("../R/weibullRound.R")
source("../R/priors.R")
source("../R/mcmcCodeInfection.R")
source("../R/mcmcRunInfection.R")

params = pandemicParams()

betaPrior=c(
 global = c(shape = 10,scale=20),
 community=c(shape=10,scale=20))

betaParams=c(global=0.5,community=0.1)

prior=pandemicPriors()

data=simEpidemicInfection(params)

output2=mcmcInfection(data,params,betaParams,prior,betaPrior,pop=10000,sigma=0.1,sigma2=0.05,runs=10)
  
 