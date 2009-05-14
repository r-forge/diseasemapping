# source in all the R files
Rfiles=system("CMD /K dir /B ..\\R\\*.R", intern=T)
Rfiles = grep("\\.R$", Rfiles, value=T)
for(D in Rfiles) {
	source(paste("..\\R\\", D, sep=""))
}	
	
# get parameters
# use the default except for hospital admission times for deadly infections
params = pandemicParams( 
  MedHospD = c(mean = 1.5, shape = 2, zeros = 0.2)
  )

# get priors
# start with default prior distributions

prior=pandemicPriors()

# set infection parameters and priors
betaParams=c(global=0.5,community=0.1)

betaPrior=c(
 global = c(shape = 10,scale=20),
 community=c(shape=10,scale=20))

# Simulate data

data=simEpidemicInfection(params,betaParams)

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
paramSampleInfection = mcmcInfection(data,params,betaParams,prior,betaPrior,pop=10000,sigma=0.1,sigma2=0.05,runs=10)
# 0.1 - Standard deviation for random walk Metropolis algorithm
# 100 - Number of iterations


# summary statistics

# get posterior mean of parameter values, transform to parameter list
paramPosteriorMean = vecParamsToList(apply(paramSample, 1, mean))

# generate hospital times
hospSample = simHospitals(paramsPosteriorSample, 100, 10)

# graph hospital times
plotHospitals(hospSample)

