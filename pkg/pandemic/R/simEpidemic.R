simEpidemic <- function(params, delta=5, days=20,
  probOnsetMissing=0.7, randomInfections = TRUE) {
  
  if(!any(names(params[["OnsMedM"]])=="scale"))
    warning("can't find scale parameter, use addScaleParameters")
  if(randomInfections) {
    infectionsPerDay=rpois(days,delta)
  } else {
    infectionsPerDay=rep(delta,days)
  }

  N=sum(infectionsPerDay)

  colNames <- c("infect", "onset",  "med", "hospital",
     "removed", "censor", "type","observedType")
  result <- matrix(NA, N, length(colNames), 
    dimnames = list(NULL, colNames))
  result <- as.data.frame(result)
  
  j=0
for(i in 1:days)
{
if(infectionsPerDay[i]>0)
{
result$infect[(j+1):(j+infectionsPerDay[i])]=i
j=j+infectionsPerDay[i]
}
}

# if all 3 probabilities provided dont depend on age
# assign states. 
  theTypes = c("M","S","D")

if(all(theTypes %in% names(params$probs))) {
  result$type = sample(factor(theTypes, levels=theTypes),
    N, replace=T, prob=params$probs)
 } else { # probs vary with age  

# construct probabilities of being S or D, which might depend on age
result$age = sample(1:80, dim(result)[1], replace=T)
result$type=NA

# create matrix to store probabilities
probMat = data.frame(probS = result$type, probD=result$type)

for(Dprob in c("S","D")) {
  if(any(names(params$ageProbs)==Dprob))  {
  # if probability corresponding to Dprob changes with age
  # assign the approrpriate probability for the age group.
  probMat[,paste("prob",Dprob,sep="")] = 
    approx(params$ageProbs[[Dprob]]$age, params$ageProbs[[Dprob]]$prob, result$age)$y
  } else {
      # if it doesnt change with age, assign the probability provided to all observations
      probMat[,paste("prob",Dprob,sep="")] = params$probs[Dprob]
  }
}
# simluate deaths
# note that now probS is prob of S conditional on not D
result$type = rbinom(dim(result)[1], 1, probMat$probD)
# convert 1's to D, 0's to NA
result$type = c(NA, "D")[result$type + 1]
#assign those which aren.t D as M or S
theNA = is.na(result$type)
theS = rbinom(sum(theNA), 1, probMat[theNA, "probS"])
result[theNA,"type"] = c("M","S")[theS+1]

} #end if probs vary with age
 

  
  
  theTypes = c(theTypes, "med","hosp")
  result$observedType = factor(result$type, levels=theTypes)
  
  result$onset=result[,"infect"]+round(rweibull(N,
    shape= getVecParams(params, "InfOns", "shape"),
    scale= getVecParams(params, "InfOns", "scale")
    ) )   
  
  result$med=result[,"onset"]+round(rweibull(N,
    shape= getVecParams(params, "OnsMed", "shape")[as.character(result$type)],
    scale= getVecParams(params, "OnsMed", "scale")[as.character(result$type)]
    ) )   

N=sum(result$med<=days)

# get rid of observations which havent had their med event before the end of 
# the study period
result=result[result$med<=days,]

  # onset times
  haveOnset = as.logical(rbinom(N, 1, 1-probOnsetMissing))
  
  result[haveOnset==F,"onset"]  = NA
  result[,"infect"]=NA

  # mild infections, 
  # recovery date
  theMild = which(result$type=="M")
  Nmild = length(theMild)
  # lost to followup
  theLost = theMild[as.logical(rbinom(Nmild, 1, params[["MedRec"]]["lost"]))]
  result[theLost,c("type", "removed")] = NA
  result[theLost,"censor"] = days
  result[theLost,"observedType"] = "med"

  # if not lost, generate recovery times

  Nmild = sum(result$observedType=="M")
  
     result[result$observedType=="M","removed"]=result[result$observedType=="M","med"]+round(rweibull(Nmild,
    shape= getVecParams(params, "MedRec", "shape"),
    scale= getVecParams(params, "MedRec", "scale")
    ) )   
  

  for(i in 1:N)
  {
  if((result[i,"observedType"]=="M")&(result[i,"removed"]>days))
  {
  result[i,"type"]=NA
  result[i,"observedType"]="med"
  result[i,"removed"]=NA
  result[i,"censor"]=days  
  }
  } 


  # serious and deadly infection
  # hospitalization
  theSerious = result$type=="S"
  theDeadly = result$type=="D"
  theNA = is.na(result$type)
  theSerious[theNA] = theDeadly[theNA] = FALSE
  
  result[theSerious, "hospital"] = result[theSerious,"med"]+
    rweibullRound(sum(theSerious), params[["MedHospS"]])
  result[theDeadly, "hospital"] = result[theDeadly,"med"]+
    rweibullRound(sum(theDeadly), params[["MedHospD"]])

  
  censored = (result$hospital) > days
  censored[is.na(censored)] = FALSE
  result[censored,c("type", "hospital")] <- NA
  result[censored,"censor"] <- days
  result[censored,"observedType"] <- "med"

  

  # recovery or death
  inHosp = !is.na(result$hospital)
  seriousHosp = theSerious & inHosp
  deadlyHosp = theDeadly & inHosp
  result[seriousHosp,"removed"] = result[seriousHosp,"hospital"]+ 
     rweibullRound(sum(seriousHosp), params[["HospRec"]])

  result[deadlyHosp,"removed"] = result[deadlyHosp,"hospital"]+ 
     rweibullRound(sum(deadlyHosp), params[["HospDeath"]])

  censored = inHosp & ((result$removed > days) )
  censored[is.na(censored)] = FALSE
  result[censored,c("type","removed")] <- NA
  result[censored,"censor"] <- days
  result[censored,"observedType"] <- "hosp"
 
  result=result[order(result$med),]

  
  result[,"onset"]=result[,"onset"]-result[,"med"]
result[,"hospital"]=result[,"hospital"]-result[,"med"]
result[,"removed"]=result[,"removed"]-result[,"med"]

  result
  
}


