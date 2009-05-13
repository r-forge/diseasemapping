simEpidemicInfection <- function(params, beta=c(global=0.5,community=0.5), pop=100000, days=20,
  probOnsetMissing=0.7) {
  
  if(!any(names(params[["OnsMedM"]])=="scale"))
    warning("can't find scale parameter, use addScaleParameters")

  colNames <- c("infect", "onset",  "med", "hospital",
     "removed", "censor", "type","observedType")
  result <- matrix(NA, nrow=0, length(colNames))
#  result <- as.data.frame(result)
  
  sus=pop

for(i in 1:days)
{
count=0

if(nrow(result)>0)
{
count=((result[,1]<i)&(result[,3]>=i))
count=sum(count)
}

sus=pop-nrow(result)

infections=rpois(1,((beta["global"]+beta["community"]*count)*sus/pop))

if(infections>0)
{

 resultx <- matrix(NA, nrow=infections, length(colNames))

  resultx[,1]=i

  theTypes = c("M","S","D")
  typex= sample(factor(theTypes, levels=theTypes),
    infections, replace=T, prob=params$probs)
  
#  theTypes = c(theTypes, "med","hosp")
#  result$observedType = factor(result$type, levels=theTypes)
  
  resultx[,2]=resultx[,1]+round(rweibull(infections,
    shape= getVecParams(params, "InfOns", "shape"),
    scale= getVecParams(params, "InfOns", "scale")
    ) )   
  
  resultx[,3]=resultx[,2]+round(rweibull(infections,
    shape= getVecParams(params, "OnsMed", "shape")[as.character(typex)],
    scale= getVecParams(params, "OnsMed", "scale")[as.character(typex)]
    ) )   

  resultx[,7]=typex 

result=rbind(result,resultx)

}

}

#result

#}

#simEpidemic(params)

N=nrow(result)

  colNames <- c("infect", "onset",  "med", "hospital",
     "removed", "censor", "type","observedType")
  result <- matrix(result, nrow=N,byrow=F, length(colNames),
dimnames = list(NULL, colNames))
result <- as.data.frame(result)
 


typeN=result[,"type"]
for(i in 1:nrow(result)) 
{
result[typeN==1,"type"]="M"
result[typeN==2,"type"]="S"
result[typeN==3,"type"]="D"
}

#result

#}

#simEpidemic(params)

N=sum(result$med<=days)


result=result[result$med<=days,]


#result

#}

#simEpidemic(params)

N=sum(result$med<=days)

  theTypes = c(theTypes, "med","hosp")
  result$observedType = factor(result$type, levels=theTypes)


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
  

#result

#}

#simEpidemic(params)


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


