paramUpdateInfection=function(params,prior,data,UInfX,name,x,sigma)
{
paramsNew=params
paramsNew[[name]][x]=abs(rnorm(1,paramsNew[[name]][x],sigma))
paramsNew=addMeanParameters(paramsNew)
ratio=Like1(data,params,paramsNew,name)

ratio=ratio*LikeUno(UInfX,params,paramsNew)

ratio=ratio*
   dprior(paramsNew[[name]]["mean"],prior[[name]]["mean"]$mean)/
  dprior(params[[name]]["mean"],prior[[name]]["mean"]$mean)
ratio=ratio*
  dprior(paramsNew[[name]]["shape"],prior[[name]]["shape"]$shape)/
  dprior(params[[name]]["shape"],prior[[name]]["shape"]$shape)
ratio= ratio*
  dprior(paramsNew[[name]]["zeros"],prior[[name]]["zeros"]$zeros)/
  dprior(params[[name]]["zeros"],prior[[name]]["zeros"]$zeros)
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

#
# Below required for mcmc3
#

InfectionNO=function(data,params)
{
Infection=data[,"med"]+data[,"infect"]
timeLag=max(data[is.na(data[,"censor"])==F,"censor"])-min(Infection)
NumberInfected=matrix(NA,ncol=3,nrow=(timeLag))
for(i in 0:timeLag)
{
NumberInfected[i,]=c(sum(((data[,"infect"]+data[,"med"])<(i+min(Infection)))&(data[,"med"]>=(i+min(Infection)))),
sum((data[,"infect"]+data[,"med"])==(i+min(Infection))),
i+min(Infection))
}
NumberInfected
}
# Column 1 -- Number of observed infectives
# Column 2 -- Number of observed infections (on day)
# Column 3 -- Day
#InfectionNO(data,params)




InfectionUnobserved=function(data,params,betaParams,pop)
{
Infection=data[,"med"]+data[,"infect"]
timeLag=max(data[is.na(data[,"censor"])==F,"censor"])-min(Infection)
sus=pop
INFE=InfectionNO(data,params)
UInfected=matrix(0,ncol=4,nrow=(timeLag))
UInfected[1,]=c(0,INFE[1,2],INFE[1,2],1+min(Infection))
for(i in 2:(timeLag))
{
UIn=0
if(i >1) UIn=sum(UInfected[1:(i-1),1])
# Need to include sus/pop
UInfected[i,]=c( rpois(1,(betaParams["global"]+betaParams["community"]*(UIn+INFE[i,1]))*
((sus-sum(UInfected[1:(i-1),2]))/pop)*
(1-(sumWeibull((timeLag-i),params$InfOns,params$OnsMedM)*params$probs["M"]+
sumWeibull((timeLag-i),params$InfOns,params$OnsMedS)*params$probs["S"]+
sumWeibull((timeLag-i),params$InfOns,params$OnsMedD)*params$probs["D"]))),
0,0,
              i+min(Infection))
UInfected[i,2]=INFE[i,2]+UInfected[i,1]
UInfected[i,3]=INFE[i,1]+sum(UInfected[(1:i),1])
}
UInfected
}
# Column 1 -- Number of unobserved infections (on day)
# Column 2 -- Total number of infections (on day)
# Column 3 -- Total number of infectives
# Column 4 -- Day


LikeUno=function(UInfX,params,paramsNew)
{
ratio=1
timeLag=UInfX[nrow(UInfX),4]
for(i in 1:nrow(UInfX))
{
bottom=1-(sumWeibull((timeLag-i),params$InfOns,params$OnsMedM)*params$probs["M"]+
sumWeibull((timeLag-i),params$InfOns,params$OnsMedS)*params$probs["S"]+
sumWeibull((timeLag-i),params$InfOns,params$OnsMedD)*params$probs["D"])
top=1-(sumWeibull((timeLag-i),paramsNew$InfOns,paramsNew$OnsMedM)*paramsNew$probs["M"]+
sumWeibull((timeLag-i),paramsNew$InfOns,paramsNew$OnsMedS)*paramsNew$probs["S"]+
sumWeibull((timeLag-i),paramsNew$InfOns,paramsNew$OnsMedD)*paramsNew$probs["D"])
ratio=ratio*(top/bottom)^UInfX[i,2]
}
ratio
}


