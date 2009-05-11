paramUpdate=function(params,data,name,x,sigma)
{
paramsNew=params
paramsNew[[name]][x]=abs(rnorm(1,paramsNew[[name]][x],sigma))
paramsNew=addMeanParameters(paramsNew)
ratio=Like1(data,params,paramsNew,name)
ratio=ratio*
   dprior(paramsNew[[name]["mean"],prior[[name]]["mean"]$mean)/
  dprior(params[[name]["mean"],prior[[name]]["mean"]$mean)
ratio=ratio*
  dprior(paramsNew[[name]["shape"],prior[[name]]["shape"]$shape)/
  dprior(params[[name]["shape"],prior[[name]]["shape"]$shape)
ratio= ratio*
  dprior(paramsNew[[name]["zeros"],prior[[name]]["zeros"]$zeros)/
  dprior(params[[name]["zeros"],prior[[name]]["zeros"]$zeros)
if(ratio>runif(1)) params=paramsNew
params
}

probsUpdate=function(data,probs,prior)
{
probs["D"]=rbeta(1,sum(data$type=="D")+prior$probs$fatality["shape1"],
    sum(data$type=="M")+sum(data$type=="S")+prior$probs$fatality["shape2"])

probs["S"]=(1-probs["D"])*rbeta(1,sum(data$type=="S")+prior$probs$hosp["shape1"],
    sum(data$type=="M")+prior$probs$hosp["shape2"])

probs["M"]=1-probs["D"]-probs["S"]

probs
}


Likelihood=function(data,dataNew,params,paramsNew)
{
like=1
delta=8
like=prod(dweibullzero(dataNew[,"onset"]-dataNew[,"infect"],paramsNew$InfOns))/
prod(dweibullzero(data[,"onset"]-data[,"infect"],params$InfOns))
like=like*prod(dweibullzero(-dataNew[dataNew$type=="M","onset"],paramsNew$OnsMedM))/
prod(dweibullzero(-data[data$type=="M","onset"],params$OnsMedM))
like=like*prod(dweibullzero(-dataNew[dataNew$type=="S","onset"],paramsNew$OnsMedS))/
prod(dweibullzero(-data[data$type=="S","onset"],params$OnsMedS))
like=like*prod(dweibullzero(-dataNew[dataNew$type=="D","onset"],paramsNew$OnsMedD))/
prod(dweibullzero(-data[data$type=="D","onset"],params$OnsMedD))
like=like*prod(dweibullzero(dataNew[dataNew$type=="M","removed"],paramsNew$MedRec))/
prod(dweibullzero(data[data$type=="M","removed"],params$MedRec))
like=like*prod(dweibullzero(dataNew[dataNew$type=="S","hospital"],paramsNew$MedHospS))/
prod(dweibullzero(data[data$type=="S","hospital"],params$MedHospS))
like=like*prod(dweibullzero(dataNew[dataNew$type=="D","hospital"],paramsNew$MedHospD))/
prod(dweibullzero(data[data$type=="D","hospital"],params$MedHospD))
like=like*prod(dweibullzero(dataNew[dataNew$type=="S","removed"],paramsNew$HospRec))/
prod(dweibullzero(data[data$type=="S","removed"],params$HospRec))
like=like*prod(dweibullzero(dataNew[dataNew$type=="D","removed"],paramsNew$HospDeath))/
prod(dweibullzero(data[data$type=="D","removed"],params$HospDeath))

like=like*paramsNew$probs["M"]^sum(dataNew$type=="M")*
paramsNew$probs["S"]^sum(dataNew$type=="S")*
paramsNew$probs["D"]^sum(dataNew$type=="D")/
(params$probs["M"]^sum(data$type=="M")*
params$probs["S"]^sum(data$type=="S")*
params$probs["D"]^sum(data$type=="D"))

Eprime=InfectionP(dataNew,paramsNew)
E=InfectionP(data,params)
like=like*(Eprime/E)^nrow(data)*exp(-delta*(Eprime-E))

like
}


Like1=function(data,params,paramsNew,name)
{
like=1
if(name=="InfOns")
like=prod(dweibullzero(dataNew[,"onset"]-dataNew[,"infect"],paramsNew$InfOns))/
prod(dweibullzero(data[,"onset"]-data[,"infect"],params$InfOns))
if(name=="OnsMedM")
like=like*prod(dweibullzero(-dataNew[dataNew$type=="M","onset"],paramsNew$OnsMedM))/
prod(dweibullzero(-data[data$type=="M","onset"],params$OnsMedM))
if(name=="OnsMedS")
like=like*prod(dweibullzero(-dataNew[dataNew$type=="S","onset"],paramsNew$OnsMedS))/
prod(dweibullzero(-data[data$type=="S","onset"],params$OnsMedS))
if(name=="OnsMedD")
like=like*prod(dweibullzero(-dataNew[dataNew$type=="D","onset"],paramsNew$OnsMedD))/
prod(dweibullzero(-data[data$type=="D","onset"],params$OnsMedD))
if(name=="MedRec")
like=like*prod(dweibullzero(dataNew[dataNew$type=="M","removed"],paramsNew$MedRec))/
prod(dweibullzero(data[data$type=="M","removed"],params$MedRec))
if(name=="MedHospS")
like=like*prod(dweibullzero(dataNew[dataNew$type=="S","hospital"],paramsNew$MedHospS))/
prod(dweibullzero(data[data$type=="S","hospital"],params$MedHospS))
if(name=="MedHospD")
like=like*prod(dweibullzero(dataNew[dataNew$type=="D","hospital"],paramsNew$MedHospD))/
prod(dweibullzero(data[data$type=="D","hospital"],params$MedHospD))
if(name=="HospRec") 
like=like*prod(dweibullzero(dataNew[dataNew$type=="S","removed"],paramsNew$HospRec))/
prod(dweibullzero(data[data$type=="S","removed"],params$HospRec))
if(name=="HospDeath")
like=like*prod(dweibullzero(dataNew[dataNew$type=="D","removed"],paramsNew$HospDeath))/
prod(dweibullzero(data[data$type=="D","removed"],params$HospDeath))
like
}



InfectionP=function(data,params)
{
Infection=data[,"med"]+data[,"infect"]
Removed=data[,"med"]+data[,"removed"]
Removed[is.na(data[,"hospital"])==FALSE]=Removed[is.na(data[,"hospital"])==FALSE]+
data[is.na(data[,"hospital"])==FALSE,"hospital"]
timeLag=max(data[is.na(data[,"censor"])==F,"censor"])-min(Infection)
ExpectedNoCases=0
for(i in 0:timeLag)
{
ExpectedNoCases=ExpectedNoCases+sumWeibull(i,params$InfOns,params$OnsMedM)*params$probs["M"]+
sumWeibull(i,params$InfOns,params$OnsMedS)*params$probs["S"]+
sumWeibull(i,params$InfOns,params$OnsMedD)*params$probs["D"]
}
ExpectedNoCases
}


Like2=function(data,dataNew,params,i)
{
like=1
delta=8
like=dweibullzero(dataNew[i,"onset"]-dataNew[i,"infect"],params$InfOns)/
dweibullzero(data[i,"onset"]-data[i,"infect"],params$InfOns)
like=like*Like3(dataNew,i)/Like3(data,i)

like=like*params$probs[dataNew[i,"type"]]/params$probs[data[i,"type"]]

Eprime=InfectionP(dataNew,params)
E=InfectionP(data,params)
like=like*(Eprime/E)^nrow(data)*exp(-delta*(Eprime-E))

like
}

Like3=function(data,i)
{
like=1
name=data[i,"type"]
if(name=="M")
{
like=like*dweibullzero(-data[i,"onset"],params$OnsMedM)
like=like*dweibullzero(data[i,"removed"],params$MedRec)
}
if(name=="S")
{
like=like*dweibullzero(-data[i,"onset"],params$OnsMedS)
like=like*dweibullzero(data[i,"hospital"],params$MedHospS)
like=like*dweibullzero(data[i,"removed"],params$HospRec)
}
if(name=="D")
{
like=like*dweibullzero(-data[i,"onset"],params$OnsMedD)
like=like*dweibullzero(data[i,"hospital"],params$MedHospD)
like=like*dweibullzero(data[i,"removed"],params$HospDeath)
}
like
}

priorEval=function(params,prior,name)
{
prob=dprior(2.5,prior$name$shape)