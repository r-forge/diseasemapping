dataAugment = function(data, params) {

# create probD and probS columns
for(Dprob in c("D","S") ) {
  probCol = paste("prob", Dprob, sep="")
  if(Dprob %in% names(params$probs)) {
    data[,probCol] = params$probs[Dprob] 
  } else {
    data[,probCol] = approx(params$ageProbs[[Dprob]]$age, 
      params$ageProbs[[Dprob]]$prob, data[,"age"])$y  # for age varying probabilities
  }
}
#  convert 2 conditional probabilities                                             # into 3 marginal probabilities
if(!all(c("M","S","D") %in% names(params$probs))) {        # do this if given conditional probabilities of S and D, may need this because if probabilities vary with age (ageProbs), will only need S and D probs
  data$probS = (1-data$probD)*data$probS   #marginal prob of being in serious state = (1-probD)*probS, where probD = marginal probabiliity of being in deadly state, probS = P(hospitalization|not dying)
  data$probM = 1 - (data$probD + data$probS)   # marginal prob of being in mild state = 1 - (prob of being in deadly state + prob of being in serious state)
} else {                                  # do this if we are given marginal probabilities of M, S, D
  data$probM = params$probs["M"]  # define probM to be the marginal probability of being in the mild state
}
    # probD = marginal prob of dying
    # probS = p(hospitalization|not dying)
    # probM = mild infection type

# hosps
inHosp = which(data$observedType == "hosp")
if(length(inHosp)) {
inHospTimes = data[inHosp, "censor"] - data[inHosp, "hospital"] - data[inHosp, "med"]
probCensorGivenDeadly = 1-pweibullRound(inHospTimes, params$HospDeath)
probCensorGivenSerious = 1-pweibullRound(inHospTimes, params$HospRec)

probCensorAndDeadly = probCensorGivenDeadly * data[inHosp,"probD"] 

probDeadlyGivenCensor = 
probCensorAndDeadly / 
  (probCensorAndDeadly + probCensorGivenSerious*data[inHosp,"probS"]) 
# this is all conditional on being in hospital.
# should divide each term by pr(in hospital) but they all cancel out.
  
  
 if(any(is.na(probDeadlyGivenCensor))) {
  warning("problem with data augmentation, saving data and params in debugstuff.RData")
  save(params, data, file="debugstuff.RData")
 } 
  

  
#if(length(inHosp)) { 
data[inHosp,"type"] = c("S","D")[1+rbinom(length(inHosp), 1, probDeadlyGivenCensor)]  
}

# meds
inMed = which(data$observedType == "med")
if(length(inMed)) {
inMedTimes = data[inMed, "censor"] - data[inMed, "med"]

# note that these probabilities are only proportional
# havent divided by pr(being med at time X)
# haven't divided by the probabililty given the time we're at (the data)
probOf = data.frame(
  "L" = data[inMed,"probM"] * params$MedRec["lost"], # prob of mild and lost | medical at time x
   M= (1-pweibullRound(inMedTimes, params$MedRec)) *
      (data[inMed,"probM"]* (1- params$MedRec["lost"]) ) , # prob of mild and not lost | medical at time x 
  "S"= data[inMed,"probS"]*
    (1-pweibullRound(inMedTimes, params$MedHospS)), # prob of serious | medical at time x
  "D"= data[inMed,"probD"]*
    (1-pweibullRound(inMedTimes, params$MedHospD) ) # prob of deadly | medical at time x
  )
 
# generate states
if(any(is.na(probOf))) {
  cat("NA probabilities resulting from these parameters")
  print(params)
  }

probOf[is.na(probOf)] = 0

states = apply(probOf, 1, function(qq) {   
    sample(colnames(probOf), 1, prob=qq)
} )

theLost = states=="L" 
states[theLost]="M"
data$lost = FALSE
data[inMed[theLost],"lost"] = T
data[inMed,"type"] = states

} # end if length(inMed)

 #generate times of medical visits
 
 needMed = is.na(data$med)
if(any(needMed)) { 
 if(any(is.na(data[needMed,"infect"])) ) {
    warning("some subjects where infection times and medical visit time are both missing")
 }

 # simulate onset times   
  needOnset=is.na(data$onset)& needMed  
  if(any(needOnset)) {

  data[needOnset,"onset"]  = data[needOnset, "infect"] + 
    round(rweibull(sum(needOnset),
        shape= getVecParams(params, "InfOns", "shape"),
        scale= getVecParams(params, "InfOns", "scale")
      ) )   

    }   # end if any needOnset


   needMedCensor = needMed & !is.na(data$censor)
   for(Dtype in c("M","S","D")) {
    needHere = needMedCensor & data$type==Dtype
    if(any(needHere)) {
  data[needHere,"med"] =    data[needHere, "onset"] + 
   censorweibull(           # get caught in an infinite loop!
    getVecParams(params, paste("OnsMed",Dtype, sep=""), "shape"),
    getVecParams(params, paste("OnsMed",Dtype, sep=""), "scale"),
(data[needHere,"censor"]-data[needHere,"onset"])
    ) 
    
    }
    
    }
   
  needMedNotC = needMed & is.na(data$censor)
  
  data[needMedNotC,"med"]  = data[needMedNotC, "onset"] + 
    round(rweibull(sum(needMedNotC),
      shape= getVecParams(params, "OnsMed", "shape")[as.character(data[needMedNotC,"type"])],
    scale= getVecParams(params, "OnsMed", "scale")[as.character(data[needMedNotC,"type"])]
     ) )   
   
   data[needMed,"onset"] =      data[needMed,"onset"]-data[needMed,"med"] 
   data[needMed,"infect"] =       data[needMed,"infect"] - data[needMed,"med"]  -   data[needMed,"onset"]
   




 } # end if length needMed
 
  needOnset=is.na(data$onset)
  if(any(needOnset)) {
  data[needOnset,"onset"]  = - round(rweibull(sum(needOnset),
    shape= getVecParams(params, "OnsMed", "shape")[as.character(data$type)],
    scale= getVecParams(params, "OnsMed", "scale")[as.character(data$type)]
    ) )   
    }

  needInfect=is.na(data$infect)
  if(any(needInfect)) {
  data[needInfect,"infect"]  = data[needInfect,"onset"]- round(rweibull(sum(needInfect),
    shape= getVecParams(params, "InfOns", "shape"),
    scale= getVecParams(params, "InfOns", "scale")
    ) ) 
    }

 
# serious, need hospital, censored 

   needhospital = (data$type =="S") & is.na(data$hospital) &(!is.na(data$censor))
  # but for censoring timme to be relevant, medical visit is before censoring time
   needhospital = needhospital & (data$med < data$censor)
   needhospital[is.na(needhospital)] = F

if (any(needhospital))
{
   data[needhospital,"hospital"]=censorweibull(
    getVecParams(params, "MedHospS", "shape"),
    getVecParams(params, "MedHospS", "scale"),
(data[needhospital,"censor"]-data[needhospital,"med"])
    ) 
}   
   # need hospital, serious, not censored.
   needhospital = (data$type =="S") & is.na(data$hospital) & (!needhospital)
   
    if (any(needhospital))
{
   data[needhospital,"hospital"]=rweibullRound( sum(needhospital),
  params$MedHospS)
     
   

}
          # deadly, need hospotal, censored
   needhospital = (data$type =="D") & is.na(data$hospital) &(!is.na(data$censor))
     needhospital = needhospital & (data$med < data$censor)
     needhospital[is.na(needhospital)] = F

if (sum(needhospital)>0)
{
   data[needhospital,"hospital"]=censorweibull(
     getVecParams(params, "MedHospD", "shape"),
    getVecParams(params, "MedHospD", "scale"),
(data[needhospital,"censor"]-data[needhospital,"med"])
    ) 

   
}
   # need hospital, deadly, not censored.
   needhospital = (data$type =="D") & is.na(data$hospital) &(!needhospital)
    if (any(needhospital)) {
   data[needhospital,"hospital"]=rweibullRound( sum(needhospital),
    params$MedHospD)
}



needremoved = data$type=="M" & is.na(data$removed)
   # depends if you're lost or not


#if not lost but censored
    needremovedC = !data[,"lost"]  & !is.na(data$censor)  & needremoved 


if (sum(needremovedC)>0)
     data[needremovedC,"removed"]=censorweibull(
getVecParams(params, "MedRec", "shape"),
 getVecParams(params, "MedRec", "scale"),
(data[needremovedC,"censor"]-data[needremovedC,"med"])
    )   


#lost or uncensored and need removed
lostNeedRemoved =  needremoved & ! needremovedC

if (any(lostNeedRemoved) )
   data[lostNeedRemoved,"removed"]=round(rweibull(sum(lostNeedRemoved),
  shape= getVecParams(params, "MedRec", "shape"),
    scale= getVecParams(params, "MedRec", "scale")
    ) )   
    


    # recovery from hospital , censoring
    
   needremoved= data$type=="S" & is.na(data$removed) & (!is.na(data$censor))

if (sum(needremoved)>0)
   data[needremoved,"removed"]=censorweibull(
getVecParams(params, "HospRec", "shape"),
 getVecParams(params, "HospRec", "scale"),
(data[needremoved,"censor"]-data[needremoved,"hospital"]-data[needremoved,"med"])
    ) 

    # same but not censored
   needremoved= data$type=="S" & is.na(data$removed) & is.na(data$censor)

if (sum(needremoved)>0)
   data[needremoved,"removed"]=rweibullRound(sum(needremoved),
    params$HospRec)



   # death dates, censoring
   needremoved= data$type=="D" & is.na(data$removed) & (!is.na(data$censor))

if (sum(needremoved)>0)
   data[needremoved,"removed"]=censorweibull(
getVecParams(params, "HospDeath", "shape"),
 getVecParams(params, "HospDeath", "scale"),
(data[needremoved,"censor"]-data[needremoved,"hospital"]-data[needremoved,"med"])
    ) 

    # same, not censored
     needremoved= data$type=="D" & is.na(data$removed) & is.na(data$censor)

if (sum(needremoved)>0)
   data[needremoved,"removed"]=rweibullRound(sum(needremoved),
    params$HospDeath)


data

}

