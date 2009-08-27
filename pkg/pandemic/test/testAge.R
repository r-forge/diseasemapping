library(DPpackage)
library(pandemic)


prior= pandemicPriors(probs=psProbPriors())


writePrior(prior, file="testprior.txt")
prior2 =readPrior("testprior.txt")


params = pandemicParams(
  probs=NULL,
  ageProbs=list(
    S=addAgeProbs(),
    D=addAgeProbs()
  )
)

x = simEpidemic(params)

MCMCscaling = mcmcScale(params, sigma=0.5, minScale = 0.01, maxScale = 0.2)
MCMCscaling$HospDeath["mean"] = 0.1

postSample = mcmcPandemic(x,params, prior, sigma=MCMCscaling, runs=20, thin=10)

plotAgeProbs(postSample, type="fatality", 
  quantiles=c(0.1, 0.9),prior=NULL) 


theprior = list(
beta0=log(priorProb/(1-priorProb)), # prior mean on the logit scale 
#Sbeta0=.5^2, # variance matrix
  taub1=2, taub2=.1/2) # inverse of penalty parameters, so small means smooth    

 theprior$Sbeta0 = ( (upperLogit - theprior$beta0)/2  )

   
   deathfit <-PSgam(formula=data$death~ps(data$age,k=Nsplines,degree=degree,pord=pord),
                family=binomial(logit),prior=theprior,
                mcmc=PSmcmc,ngrid=Ngrid,
                state=attributes(probs$D)$state,
                status=is.null(attributes(probs$D)$state) )

plot(deathfit)


