library(pandemic)
	
# get parameters
# use the default except for hospital admission times for deadly infections
params = pandemicParams( 
  MedHospD = c(mean = 1.5, shape = 2, zeros = 0.2)
  )

# simulate data
# delta is mean number of infections per day
data = simEpidemic(params, delta=3, days=20)

# get priors
# start with default prior distributions
priors = pandemicPriors()
# write hyperparameters to a file
writePrior(priors, "hyperparameters.txt")
# make some graphs to fine better priors
plotPrior(
x= meanShapeZerosPrior(
 mean=c(mean=5, sd=0.5), 
 shape=c(mean=1.5, sd=0.2), 
 zeros = c(mean=0.2, sd=0.02))
 )
# edit the file hyperparameters.txt to change the prior distributions
#   and read it in
priors = readPrior("hyperparameters.txt")

# fit the model
paramSample = mcmcPandemic(data,params,priors,0.1,100)

# 0.1 - Standard deviation for random walk Metropolis algorithm
# 100 - Number of iterations

# summary statistics
postSummary=t(apply(paramSample, 2, function(qq) 
  c(mean=mean(qq), quantile(qq, prob=c(0.025, 0.5, 0.975)))))
postSummary = cbind(true=unlist(params)[rownames(postSummary)], postSummary)  


# graph of prior and posterior histogram
plotPrior(priors, paramSample, transition="HospDeath")


# get posterior mean of parameter values, transform to parameter list
paramPosteriorMean = vecParamsToList(postSummary[,"mean"])

# generate hospital times
hospSample = simHospitals(paramSample , 100, 10)

# graph hospital times
plotHospitals(hospSample)

