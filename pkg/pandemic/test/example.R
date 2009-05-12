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

# simulate data
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
paramSample = NA

# summary statistics

# get posterior mean of parameter values, transform to parameter list
paramPosteriorMean = vecParamsToList(apply(paramSample, 1, mean))

# generate hospital times
hospSample = simHospitals(paramsPosteriorSample, 100, 10)

# graph hospital times
plotHospitals(hospSample)

