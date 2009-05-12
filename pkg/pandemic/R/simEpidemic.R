simEpidemic <- function(params, delta=5, days=20,
  probOnsetMissing=0.7) {
  
  if(!any(names(params[["OnsMedM"]])=="scale"))
    warning("can't find scale parameter, use addScaleParameters")
  
  infectionsPerDay=rpois(days,delta)

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
#  result

#}

#simEpidemic(params)


  theTypes = c("M","S","D")
  result$type = sample(factor(theTypes, levels=theTypes),
    N, replace=T, prob=params$probs)
  
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

  result
  
}


