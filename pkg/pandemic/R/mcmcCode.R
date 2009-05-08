paramUpdate=function(param,x,sigma,periods)
{
paramNew=param
paramNew[x]=abs(rnorm(1,param[x],sigma))
ratio=prod(dweibullRound(periods,paramNew))/prod(dweibullRound(periods,param))
if(ratio>runif(1)) param=paramNew
param
}





Likelihood=function(data,dataNew,params,paramsNew)
{
like=1
delta=8
like=prod(dweibullRound(dataNew[,"onset"]-dataNew[,"infect"],paramsNew$InfOns))/
prod(dweibullRound(data[,"onset"]-data[,"infect"],params$InfOns))
like=like*prod(dweibullRound(-dataNew[dataNew$type=="M","onset"],paramsNew$OnsMedM))/
prod(dweibullRound(-data[data$type=="M","onset"],params$OnsMedM))
like=like*prod(dweibullRound(-dataNew[dataNew$type=="S","onset"],paramsNew$OnsMedS))/
prod(dweibullRound(-data[data$type=="S","onset"],params$OnsMedS))
like=like*prod(dweibullRound(-dataNew[dataNew$type=="D","onset"],paramsNew$OnsMedD))/
prod(dweibullRound(-data[data$type=="D","onset"],params$OnsMedD))
like=like*prod(dweibullRound(dataNew[dataNew$type=="M","removed"],paramsNew$MedRec))/
prod(dweibullRound(data[data$type=="M","removed"],params$MedRec))
like=like*prod(dweibullRound(dataNew[dataNew$type=="S","hospital"],paramsNew$MedHospS))/
prod(dweibullRound(data[data$type=="S","hospital"],params$MedHospS))
like=like*prod(dweibullRound(dataNew[dataNew$type=="D","hospital"],paramsNew$MedHospD))/
prod(dweibullRound(data[data$type=="D","hospital"],params$MedHospD))
like=like*prod(dweibullRound(dataNew[dataNew$type=="S","removed"],paramsNew$HospRec))/
prod(dweibullRound(data[data$type=="S","removed"],params$HospRec))
like=like*prod(dweibullRound(dataNew[dataNew$type=="D","removed"],paramsNew$HospDeath))/
prod(dweibullRound(data[data$type=="D","removed"],params$HospDeath))

like=like*paramsNew$probs["M"]^sum(dataNew$type=="M")*
paramsNew$probs["S"]^sum(dataNew$type=="S")*
paramsNew$probs["D"]^sum(dataNew$type=="D")/
(params$probs["M"]^sum(data$type=="M")*
params$probs["S"]^sum(data$type=="S")*
params$probs["D"]^sum(data$type=="D"))

Eprime=InfectionP(dataNew,paramsNew)
E=InfectionP(data,params)
like=(Eprime/E)^nrow(data)*exp(-delta*(Eprime-E))

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

InfectionP(data,params)

Likelihood(data,dataNew,params,paramsNew)



#delta=rgamma(1,(nrow(data)+1),rate=ExpectedNoCases)




