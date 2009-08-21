paramUpdate=function(params,prior,data,name,x,sigma)
{
paramsNew=params
paramsNew[[name]][x]=abs(rnorm(1,paramsNew[[name]][x],sigma[[name]][x]))
paramsNew=addMeanParameters(paramsNew)
ratio1=Like1(data,params,paramsNew,name)
ratio2=ratio1*
   dprior(paramsNew[[name]]["mean"],prior[[name]]["mean"]$mean)/
  dprior(params[[name]]["mean"],prior[[name]]["mean"]$mean)
ratio3=ratio2*
  dprior(paramsNew[[name]]["shape"],prior[[name]]["shape"]$shape)/
  dprior(params[[name]]["shape"],prior[[name]]["shape"]$shape)
ratio= ratio3*
  dprior(paramsNew[[name]]["zeros"],prior[[name]]["zeros"]$zeros)/
  dprior(params[[name]]["zeros"],prior[[name]]["zeros"]$zeros)
if(is.na(ratio)){
 cat(ratio1, ratio2, ratio3, ratio, "\n")
}

if(ratio>runif(1)) params=paramsNew
params
}

probsUpdate=function(data,probs,prior)
{

if(is.vector(probs)) {

probs["D"]=rbeta(1,sum(data$type=="D")+prior$probs$fatality["shape1"],
    sum(data$type=="M")+sum(data$type=="S")+prior$probs$fatality["shape2"])

probs["S"]=(1-probs["D"])*rbeta(1,sum(data$type=="S")+prior$probs$hosp["shape1"],
    sum(data$type=="M")+prior$probs$hosp["shape2"])

probs["M"]=(1-probs["D"])*(1-probs["S"])

} else { #DP stuff
    # NEED TO DO THIS!

}

probs

}


Like1=function(data,params,paramsNew,name)
{
like=1
if(name=="InfOns")
like=prod(dweibullzero(data[,"onset"]-data[,"infect"],paramsNew$InfOns))/
prod(dweibullzero(data[,"onset"]-data[,"infect"],params$InfOns))
if(name=="OnsMedM")
like=like*prod(dweibullzero(-data[data$type=="M","onset"],paramsNew$OnsMedM))/
prod(dweibullzero(-data[data$type=="M","onset"],params$OnsMedM))
if(name=="OnsMedS")
like=like*prod(dweibullzero(-data[data$type=="S","onset"],paramsNew$OnsMedS))/
prod(dweibullzero(-data[data$type=="S","onset"],params$OnsMedS))
if(name=="OnsMedD")
like=like*prod(dweibullzero(-data[data$type=="D","onset"],paramsNew$OnsMedD))/
prod(dweibullzero(-data[data$type=="D","onset"],params$OnsMedD))
if(name=="MedRec")
like=like*prod(dweibullzero(data[data$type=="M","removed"],paramsNew$MedRec))/
prod(dweibullzero(data[data$type=="M","removed"],params$MedRec))
if(name=="MedHospS")
like=like*prod(dweibullzero(data[data$type=="S","hospital"],paramsNew$MedHospS))/
prod(dweibullzero(data[data$type=="S","hospital"],params$MedHospS))
if(name=="MedHospD")
like=like*prod(dweibullzero(data[data$type=="D","hospital"],paramsNew$MedHospD))/
prod(dweibullzero(data[data$type=="D","hospital"],params$MedHospD))
if(name=="HospRec") 
like=like*prod(dweibullzero(data[data$type=="S","removed"],paramsNew$HospRec))/
prod(dweibullzero(data[data$type=="S","removed"],params$HospRec))
if(name=="HospDeath")
like=like*prod(dweibullzero(data[data$type=="D","removed"],paramsNew$HospDeath))/
prod(dweibullzero(data[data$type=="D","removed"],params$HospDeath))
like
}

#
# Below required for mcmc2
#

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



