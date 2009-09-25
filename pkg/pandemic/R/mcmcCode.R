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
  cat("problem with ratios in paramUpdate\n")
 cat(ratio1, ratio2, ratio3, ratio, "\n")
 ratio = 0
}

if(ratio>runif(1)) params=paramsNew
params
}

probsUpdate=function(data,probs,prior)
{

if(any(names(probs)=="M") )  { # vector of D, S, M

probs["D"]=rbeta(1,sum(data$type=="D")+prior$fatality["shape1"],
    sum(data$type=="M")+sum(data$type=="S")+prior$fatality["shape2"])

probs["S"]=(1-probs["D"])*rbeta(1,sum(data$type=="S")+prior$hosp["shape1"],
    sum(data$type=="M")+prior$hosp["shape2"])

probs["M"]=(1-probs["D"])*(1-probs["S"])

} else { #DP stuff
     PSmcmc <- list(nburn=2,
                 nsave=1,
                 nskip=0,
                 ndisplay=2)

       Ngrid = 30
       Nsplines = 6
       degree=3
       pord=1
# death
ageUniqueIndex = which(!duplicated(data$age))
ageUniqueIndex = ageUniqueIndex[order(data$age[ageUniqueIndex])]
ageUnique = data$age[ageUniqueIndex]

   data$death =data$type=="D"
   
   deathfit <-PSgam(formula=data$death~ps(data$age,k=Nsplines,degree=degree,pord=pord),
                family=binomial(logit),
                prior=prior$fatality[c("taub1","taub2","beta0","Sbeta0","tau1","tau2")],
                mcmc=PSmcmc,ngrid=Ngrid,
                state=attributes(probs$D)$state,
                status=is.null(attributes(probs$D)$state) )
# get basis functions evaluated at each of the ages in the dataset
  pred = deathfit$z[ageUniqueIndex,] %*%  deathfit$state$b + deathfit$state$beta
# or get smoothed fit evaluated on the grid?
#   pred = deathfit$save.state$pssave[,1]  + deathfit$state$beta
   pred = exp(pred) / (1+exp(pred))
   
   probs$D = data.frame(age=ageUnique, prob=pred)
   attributes(probs$D)$state = deathfit$state
   
# hosp
  data =data[!data$death,]
  ageUniqueIndex = which(!duplicated(data$age))
  ageUniqueIndex = ageUniqueIndex[order(data$age[ageUniqueIndex])]
  ageUnique = data$age[ageUniqueIndex]

  data$hosp = data$type=="S"
  hospfit <-PSgam(formula=data$hosp~ps(data$age,k=Nsplines,degree=degree,pord=pord),
                family=binomial(logit),prior=prior$hosp,
                mcmc=PSmcmc,ngrid=Ngrid,
                state=attributes(probs$S)$state,
                status=is.null(attributes(probs$S)$state))

   pred = hospfit$z[ageUniqueIndex,] %*%  hospfit$state$b + hospfit$state$beta
#   pred = hospfit$save.state$pssave[,1]  + hospfit$state$beta
   pred = exp(pred) / (1+exp(pred))
   
   if(length(ageUnique) != length(probs$D$age)) {
   # some ages in the death model aren't included in the hosp model because
   # everyone in that age group is dead
   # in that case, use interpolation to find probs at these ages
   allAges = c(ageUnique, hospfit$xreal)
   
   predFromPSgam = hospfit$save.state$pssave[,1] + hospfit$coef["(Intercept)"]
   predFromPSgam = exp(predFromPSgam) / (1+exp(predFromPSgam))
   
   allpred = c(pred, predFromPSgam)
   
   probs$S = data.frame(
    age=probs$D$age, 
    prob=approx(allAges, allpred, probs$D$age)$y
    )
   } else {
    probs$S = data.frame(age=ageUnique, prob=pred)
   }

   attributes(probs$S)$state = hospfit$state
   
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



