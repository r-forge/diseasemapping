simEpidemic <- function(params, N, days=10,
  probOnsetMissing=0.7) {
  
  if(!any(names(params[["OnsMedM"]])=="scale"))
    warning("can't find scale parameter, use addScaleParameters")
  
  colNames <- c("infect", "onset",  "med", "hospital",
     "removed", "censor", "type","observedType")
  result <- matrix(NA, N, length(colNames), 
    dimnames = list(NULL, colNames))
  result <- as.data.frame(result)
  
  result$med = sample(days, N, replace=T)
  theTypes = c("M","S","D")
  result$type = sample(factor(theTypes, levels=theTypes),
    N, replace=T, prob=params$probs)
  
  theTypes = c(theTypes, "med","hosp")
  result$observedType = factor(result$type, levels=theTypes)
  
  # onset times
  haveOnset = as.logical(rbinom(N, 1, 1-probOnsetMissing))
  
  result[haveOnset,"onset"]  = - round(rweibull(sum(haveOnset),
    shape= getVecParams(params, "OnsMed", "shape")[as.character(result$type)],
    scale= getVecParams(params, "OnsMed", "scale")[as.character(result$type)]
    ) )   
  
  # mild infections, 
  # recovery date
  theMild = result$type=="M"
  result[theMild, "removed"] = 
    rweibullRound(sum(theMild), params[["MedRec"]])
    
  # censor and observedType
  censored = (result$med + result$removed) > days
  censored[is.na(censored)] = FALSE
  result[censored,"removed"] <- NA
  result[censored,"censor"] <- days
  result[censored,"observedType"] <- "med"
  
  
  # serious and deadly infection
  # hospitalization
  theSerious = result$type=="S"
  theDeadly = result$type=="D"
  result[theSerious, "hospital"] = 
    rweibullRound(sum(theSerious), params[["MedHospS"]])
  result[theDeadly, "hospital"] = 
    rweibullRound(sum(theDeadly), params[["MedHospD"]])
  
  censored = (result$med + result$hospital) > days
  censored[is.na(censored)] = FALSE
  result[censored,"hospital"] <- NA
  result[censored,"censor"] <- days
  result[censored,"observedType"] <- "med"

  # recovery or death
  inHosp = !is.na(result$hospital)
  seriousHosp = theSerious & inHosp
  deadlyHosp = theDeadly & inHosp
  result[seriousHosp,"removed"] =
     rweibullRound(sum(seriousHosp), params[["HospRec"]])

  result[deadlyHosp,"removed"] =
     rweibullRound(sum(deadlyHosp), params[["HospDeath"]])

  censored = inHosp & ((result$med + result$hospital + result$removed > days) )
  censored[is.na(censored)] = FALSE
  result[censored,"removed"] <- NA
  result[censored,"censor"] <- days
  result[censored,"observedType"] <- "hosp"

  
  result
  
}


rweibullRound = function(N, params) {
  round(rweibull(N, shape=params[["shape"]], scale=params[["scale"]]))
}

