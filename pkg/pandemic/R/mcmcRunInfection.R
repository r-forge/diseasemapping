mcmcInfection=function(xdata,params,betaParams,prior,betaprior,pop,sigma,sigma2,runs)
{
name=c("InfOns","OnsMedM","OnsMedS","OnsMedD","MedRec","MedHospS","MedHospD","HospRec","HospDeath")
name2=c("scale","shape","zeros")
name3=c("global","community")
data=dataAugment(xdata,params)

InfX=InfectionNO(data,params)
UInfX=InfectionUnobserved(data,params,betaParams,pop)

vecParams = unlist(params) 
paramsPosteriorSample = matrix(NA, ncol=length(vecParams), nrow=0, 
  dimnames = list(NULL, names(vecParams)) )

delta=matrix(NA,nrow=runs,ncol=2)

for(k in 1:runs)
{
for(i in name)
{
for(j in name2)
{
params=paramUpdateInfection(params,prior,data,UInfX,i,j,sigma)
}
}
params$probs=probsUpdate(data,params$probs,prior)

paramsPosteriorSample = rbind(paramsPosteriorSample, 
    unlist(params)[colnames(paramsPosteriorSample)])

dataNew=data
w=sample(nrow(data),1) # Could replace by 1 by a greater number if we want to update a number of rows.
dataNew[w,]=xdata[w,]
dataNew=dataAugment(dataNew,params)

if(min(dataNew[,"infect"]+dataNew[,"med"])==min(data[,"infect"]+data[,"med"]))
{
InfZ=InfectionNO(dataNew,params)
UInfZ=UInfX
UInfZ[,2]=UInfX[,2]+InfZ[,2]-InfX[,2]
UInfZ[,3]=UInfX[,3]+InfZ[,1]-InfX[,1]

MeanXT=0
MeanXTNew=0

for(t in 1:nrow(UInfX))
{
sus=pop
if(t >1) sus=pop-sum(UInfX[(1:(t-1)),2])
MeanXT[t]=(betaParams["global"]+betaParams["community"]*UInfX[t,3])*sus/pop
sus=pop
if(t >1) sus=sus-sum(UInfZ[(1:(t-1)),2])
MeanXTNew[t]=(betaParams["global"]+betaParams["community"]*UInfZ[t,3])*sus/pop
}

ratio=prod(dpois(UInfZ[,2],MeanXTNew))/prod(dpois(UInfX[,2],MeanXT))
}
else
{


InfZ=InfectionNO(dataNew,params)
UInfZ=InfectionUnobserved(dataNew,params,betaParams,pop)



MeanXT=0
MeanXTNew=0

for(t in 1:nrow(UInfX))
{
sus=pop
if(t >1) sus=pop-sum(UInfX[(1:(t-1)),2])
MeanXT[t]=(betaParams["global"]+betaParams["community"]*UInfX[t,3])*sus/pop
}
for(t in 1:nrow(UInfZ))
{
sus=pop
if(t >1) sus=sus-sum(UInfZ[(1:(t-1)),2])
MeanXTNew[t]=(betaParams["global"]+betaParams["community"]*UInfZ[t,3])*sus/pop
}


ratio=prod(dpois(UInfZ[,2],MeanXTNew))/prod(dpois(UInfX[,2],MeanXT))

}

if(ratio>runif(1))
{
data=dataNew
InfX=InfZ
UInfX=UInfZ
}

UInfX=InfectionUnobserved(data,params,betaParams,pop)

for(j in name3)
{
betaParamsNew=betaParams
betaParamsNew[j]=rnorm(1,betaParams[j],sigma2)

for(t in 1:nrow(UInfX))
{
sus=pop
if(t >1) sus=pop-sum(UInfX[(1:(t-1)),2])
MeanXT[t]=(betaParams["global"]+betaParams["community"]*UInfX[t,3])*sus/pop
MeanXTNew[t]=(betaParamsNew["global"]+betaParamsNew["community"]*UInfX[t,3])*sus/pop
}


ratio=prod(dpois(UInfX[,2],MeanXTNew))/prod(dpois(UInfX[,2],MeanXT))

ratio=ratio*dgamma(betaParamsNew["global"],betaprior["global.shape"],betaprior["global.scale"])/
dgamma(betaParams["global"],betaprior["global.shape"],betaprior["global.scale"])
ratio=ratio*dgamma(betaParamsNew["community"],betaprior["community.shape"],betaprior["community.scale"])/
dgamma(betaParams["community"],betaprior["community.shape"],betaprior["community.scale"])

if(ratio>runif(1)) betaParams=betaParamsNew
}

delta[k,]=betaParams

}
paramsPosteriorSample=cbind(paramsPosteriorSample,delta)

paramsPosteriorSample
}


