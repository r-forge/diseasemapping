mcmcScale = function(params, sigma, minScale = 0.01, maxScale = 0.2){

# dont scale probabilities which change with age.
paramsToScale = names(params)
paramsToScale = paramsToScale[!paramsToScale %in% "ageProbs"]

params=params[paramsToScale]

  for(Dtrans in paramsToScale){
    for(Dpar in names(params[[Dtrans]])){
      params[[Dtrans]][Dpar] =
        pmax(minScale,pmin(maxScale,params[[Dtrans]][Dpar] * sigma) )
    }
  }
params
}

lostUpdate = function(data, prior) {
    x = table(data[data$type=="M","lost"])
    
    rbeta(1, x["TRUE"] + prior["shape1"],x["FALSE"] + prior["shape2"])
    
}


mcmcPandemic=function(xdata,params,prior,sigma,runs, thin=1)
{
name=c("InfOns","OnsMedM","OnsMedS","OnsMedD","MedRec","MedHospS","MedHospD","HospRec","HospDeath")
name2=c("scale","shape","zeros")

if(any(names(params)=="probs")) {
    probName <- "probs"
} else {
    probName <- "ageProbs"
  # do one iteration to get ages at which probabilities are evaluated
     data=dataAugment(xdata,params)
     cat("a great many blank lines are now being printed")
     params[[probName]]=probsUpdate(data,params[[probName]],prior$probs)
} 

if(!is.list(sigma)) {
  sigma = mcmcScale(params, sigma)    
}

vecParams = unlist(params)
vecParams = vecParams[-grep("^ageProbs[[:graph:]]+age[[:digit:]]+$", names(vecParams))] 
paramsPosteriorSample = matrix(NA, ncol=length(vecParams), nrow=0, 
  dimnames = list(NULL, names(vecParams)) )

Sthin = 1:thin

for(k in 1:runs) {
  for(Dthin in Sthin) {
 
  data=dataAugment(xdata,params)
    for(i in name) {
      for(j in name2) {
        params=paramUpdate(params,prior,data,i,j,sigma)
      }
    }
  # update probabilities  
  params[[probName]]=probsUpdate(data,params[[probName]],prior$probs)

# update lost
  params$MedRec["lost"] = lostUpdate(data, prior$MedRec$lost)

  } # end Dthin
  # save parameters
  paramsPosteriorSample = rbind(paramsPosteriorSample, 
    unlist(params)[colnames(paramsPosteriorSample)])
  # save age dependent probabilities
  
}
paramsPosteriorSample
}

mcmc2=function(xdata,params,delta,prior,sigma,runs,deltaprior=c(shape=25,scale=5))
{
name=c("InfOns","OnsMedM","OnsMedS","OnsMedD","MedRec","MedHospS","MedHospD","HospRec","HospDeath")
name2=c("scale","shape","zeros")
data=dataAugment(xdata,params)

InfP=InfectionP(data,params)

vecParams = unlist(params) 
paramsPosteriorSample = matrix(NA, ncol=length(vecParams), nrow=0, 
  dimnames = list(NULL, names(vecParams)) )

for(k in 1:runs)
{
for(i in name)
{
for(j in name2)
{
params=paramUpdate(params,prior,data,i,j,sigma)
}
}
params$probs=probsUpdate(data,params$probs,prior)

delta[k]=rgamma(1,nrow(data)+deltaprior["shape"],InfP+deltaprior["scale"])


paramsPosteriorSample = rbind(paramsPosteriorSample, 
    unlist(params)[colnames(paramsPosteriorSample)])

dataNew=dataAugment(xdata,params)

InfPNew=InfectionP(dataNew,params)

ratio=(InfPNew/InfP)^nrow(data)*
  exp(-delta[k]*(InfPNew-InfP))

if(ratio>runif(1))
{
data=dataNew
InfP=InfPNew
}

}
paramsPosteriorSample=cbind(paramsPosteriorSample,delta)

paramsPosteriorSample
}

mcmc3=function(xdata,params,betaParams,prior,sigma,runs,betaprior)
{
name=c("InfOns","OnsMedM","OnsMedS","OnsMedD","MedRec","MedHospS","MedHospD","HospRec","HospDeath")
name2=c("scale","shape","zeros")
data=dataAugment(xdata,params)

InfP=InfectionP(data,params)

vecParams = unlist(params) 
paramsPosteriorSample = matrix(NA, ncol=length(vecParams), nrow=0, 
  dimnames = list(NULL, names(vecParams)) )

for(k in 1:runs)
{
for(i in name)
{
for(j in name2)
{
params=paramUpdate(params,prior,data,i,j,sigma)
}
}
params$probs=probsUpdate(data,params$probs,prior)

paramsPosteriorSample = rbind(paramsPosteriorSample, 
    unlist(params)[colnames(paramsPosteriorSample)])

dataNew=data
w=sample(nrow(data),1) # Could replace by 1 by a greater number if we want to update a number of rows.
dataNew[w,]=xdata[w,]
dataNew=dataAugment(dataNew,params)

InfPNew=InfectionP(dataNew,params)

# commented out below because delta is undefined.
ratio=(InfPNew/InfP)^nrow(data) #* exp(-delta[k]*(InfPNew-InfP))

if(ratio>runif(1))
{
data=dataNew
InfP=InfPNew
}

}
#commented out this line because delta is undefined
#paramsPosteriorSample=cbind(paramsPosteriorSample,delta)

paramsPosteriorSample
}

