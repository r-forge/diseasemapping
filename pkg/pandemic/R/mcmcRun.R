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

