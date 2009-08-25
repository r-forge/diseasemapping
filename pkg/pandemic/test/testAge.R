library(DPpackage)
library(pandemic)

  ?PSgam

priors= pandemicPriors(probs=psProbPriors())

writePrior(priors)
priors2 = readPrior("priors.txt")


params = pandemicParams(
  probs=NULL,
  ageProbs=list(
    hosp=addAgeProbs(),
    fatality=addAgeProbs()
  )
)
help files for addAgeProbs, psParams