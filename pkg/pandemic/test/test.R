source("../R/parameters.R")
source("../R/simEpidemic.R")
source("../R/weibullRound.R")
source("../R/sampleProbs.R")
source("../R/priors.R")
source("../R/readPriors.R")
  

  
params = pandemicParams()
priors = pandemicPriors()
data = simEpidemic(params, 30)
vecParams = unlist(params) 
paramsPosteriorSample = matrix(NA, ncol=length(vecParams), nrow=0, 
  dimnames = list(NULL, names(vecParams)) )



listParams = vecParamsToList(vecParams)


beds = simHospitals(paramsPosteriorSample, 100, 10)
  
source("../R/plotPriors.R")
source("../R/priors.R")
source("../R/weibullRound.R")
plotPrior(
x= meanShapeZerosPrior(
 mean=c(mean=5, sd=0.5), 
 shape=c(mean=1.5, sd=0.2), 
 zeros = c(mean=0.2, sd=0.02))
 )
  
  
  