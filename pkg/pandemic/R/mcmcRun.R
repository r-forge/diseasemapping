mcmc=function(xdata,params,prior,sigma,runs)
{
name=c("InfOns","OnsMedM","OnsMedS","OnsMedD","MedRec","MedHospS","MedHospD","HospRec","HospDeath")
name2=c("scale","shape","zeros")
data=dataAugment(xdata,params)

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

data=dataAugment(xdata,params)
}
paramsPosteriorSample
}

mcmc2=function(xdata,params,delta,prior,deltaprior=c(shape=25,scale=5),sigma,runs)
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

