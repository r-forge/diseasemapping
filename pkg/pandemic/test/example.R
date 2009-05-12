# source in all the R files
Rfiles=system("CMD /K dir /B ..\\..\\pkg\\pandemic\\R\\*.R", intern=T)
Rfiles = grep("\\.R$", Rfiles, value=T)
for(D in Rfiles) {
	source(paste("..\\..\\pkg\\pandemic\\R\\", D, sep=""))
	
# get parameters
params = pandemicParams()

# simulate data
data = simEpidemic(params, 30)

# get priors
priors = pandemicPriors()

# plot the priors
 plotPrior(priors$InfOns)

# fit the model
paramSample = NA

# summaries of parameter posteriors


# generate hospital times
hospSample = simHospitals(paramsPosteriorSample, 100, 10)

# graph hospital times
plotHospitals(hospSample)