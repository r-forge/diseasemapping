library(DPpackage)
library(pandemic)

  ?PSgam

prior= pandemicPriors(probs=psProbPriors())

writePrior(prior)
priors2 = readPrior("priors.txt")


params = pandemicParams(
  probs=NULL,
  ageProbs=list(
    S=addAgeProbs(),
    D=addAgeProbs()
  )
)

params2 = pandemicParams()


# simulate an epidemic
#delta=5; days=20;  probOnsetMissing=0.7; randomInfections = TRUE
x = simEpidemic(params)

MCMCscaling = mcmcScale(params, sigma=0.5, minScale = 0.01, maxScale = 0.2)
MCMCscaling$HospDeath["mean"] = 0.1

temp = mcmcPandemic(x,params, prior, sigma=MCMCscaling, runs=20, thin=10)

