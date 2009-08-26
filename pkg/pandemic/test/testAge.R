library(DPpackage)
library(pandemic)

  ?PSgam

priors= pandemicPriors(probs=psProbPriors())

writePrior(priors)
priors2 = readPrior("priors.txt")


params = pandemicParams(
  probs=NULL,
  ageProbs=list(
    S=addAgeProbs(),
    D=addAgeProbs()
  )
)



# simulate an epidemic
x =simEpidemic(params)

# data augmentation
xAug = dataAugment(x, params)



data=xAug
probs=params$ageProbs
prior=priors$probs

temp=probsUpdate(xAug, params$ageProbs, priors$probs)


pred = temp$save.state$pssave[,1]+ temp$state$beta
pred = exp(pred) / (1+exp(pred))
plot(temp$xreal,pred )

pred = temp$z %*%  temp$state$b + temp$state$beta
pred = exp(pred) / (1+exp(pred))


plot(xAug$age, pred )