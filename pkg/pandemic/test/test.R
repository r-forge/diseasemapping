source("../R/parameters.R")
source("../R/simEpidemic.R")
source("../R/weibullRound.R")

params = pandemicParams()
 data = simEpidemic(params, 30)
 summary(data)

 xseq=0:20

plot(xseq, dweibullRound(xseq, params$HospRec))
sum( dweibullRound(0:100, params$HospRec))

priors = pandemicPriors()
names(priors)
names(priors$InfOns)
priors$InfOns$mean
attributes(priors$InfOns$mean)$distribution

xseq = seq(0, 10, len=100)
plot(xseq, dgamma(xseq, 
  shape=priors$InfOns$mean["shape"],
  scale=priors$InfOns$mean["scale"])
  )

dprior(2, priors$InfOns$mean)
dprior(2, priors$InfOns$mean, prefix="r")

dprior(0.2, priors$InfOns$zeros)


sampleProbs(data$type, priors$probs)
  
  
  