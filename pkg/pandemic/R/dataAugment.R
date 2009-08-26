dataAugment = function(data, params) {

# create probD and probS columns
for(Dprob in c("D","S") ) {
  probCol = paste("prob", Dprob, sep="")
  if(Dprob %in% names(params$probs)) {
    data[,probCol] = params$probs["D"]
  } else {
    data[,probCol] = approx(params$ageProbs[[Dprob]]$age, 
      params$ageProbs[[Dprob]]$prob, data[,"age"])$y
  }
}
# unless three probabilities provided, convert 2 conditional probabilities
# into 3 marginal probabilities
if(!all(c("M","S","D") %in% names(params$probs))) {
  data$probS = (1-data$probD)*data$probS
  data$probM = 1 - (data$probD + data$probS)
} else {
  data$probM = params$probs["M"]
}


# hosps
inHosp = which(data$observedType == "hosp")
inHospTimes = data[inHosp, "censor"] - data[inHosp, "hospital"]
probCensorGivenDeadly = 1-pweibullRound(inHospTimes, params$HospDeath)
probCensorGivenSerious = 1-pweibullRound(inHospTimes, params$HospRec)

probDeadlyGivenCensor = probCensorGivenDeadly *data[inHosp,"probD"]

probDeadlyGivenCensor = probDeadlyGivenCensor / 
  (probDeadlyGivenCensor + probCensorGivenSerious)
  
data[inHosp,"type"] = c("S","D")[1+rbinom(length(inHosp), 1, probDeadlyGivenCensor)]  

# meds
inMed = which(data$observedType == "med")
inMedTimes = data[inMed, "censor"] - data[inMed, "med"]

# note that these probabilities are only proportional
# havent divided by pr(time)
probOf = data.frame(
  "L" = data[inMed,"probM"] * params$MedRec["lost"], # *1
   M= (1-pweibullRound(inMedTimes, params$MedRec)) *
      (data[inMed,"probM"]* (1- params$MedRec["lost"]) ) ,
  "S"= data[inMed,"probS"]*
    (1-pweibullRound(inMedTimes, params$MedHospS)),
  "D"= data[inMed,"probD"]*
    (1-pweibullRound(inMedTimes, params$MedHospD) )
  )

# generate states
states = apply(probOf, 1, function(qq) {
    sample(colnames(probOf), 1, prob=qq)
} )

theLost = states=="L" 
states[theLost]="M"
data$lost = FALSE
data[inMed[theLost],"lost"] = T
data[inMed,"type"] = states

  needOnset=is.na(data$onset)

  data[needOnset,"onset"]  = - round(rweibull(sum(needOnset),
    shape= getVecParams(params, "OnsMed", "shape")[as.character(data$type)],
    scale= getVecParams(params, "OnsMed", "scale")[as.character(data$type)]
    ) )   

  needInfect=is.na(data$infect)

  data[needInfect,"infect"]  = data[needInfect,"onset"]- round(rweibull(sum(needInfect),
    shape= getVecParams(params, "InfOns", "shape"),
    scale= getVecParams(params, "InfOns", "scale")
    ) ) 

 
   needhospital=((as.character(data$observedType)=="med")&(as.character(data$type)=="S"))

if (sum(needhospital)>0)
{
   data[needhospital,"hospital"]=censorweibull(
    getVecParams(params, "MedHospS", "shape"),
    getVecParams(params, "MedHospS", "scale"),
(data[needhospital,"censor"]-data[needhospital,"med"])
    ) 
   
  data[needhospital,"removed"]  =  round(rweibull(sum(needhospital),
  shape= getVecParams(params, "HospRec", "shape"),
    scale= getVecParams(params, "HospRec", "scale")
    ) )   
}

   needhospital=((as.character(data$observedType)=="med")&(as.character(data$type)=="D"))

if (sum(needhospital)>0)
{
   data[needhospital,"hospital"]=censorweibull(
     getVecParams(params, "MedHospD", "shape"),
    getVecParams(params, "MedHospD", "scale"),
(data[needhospital,"censor"]-data[needhospital,"med"])
    ) 

   
  data[needhospital,"removed"]  =  round(rweibull(sum(needhospital),
  shape= getVecParams(params, "HospDeath", "shape"),
    scale= getVecParams(params, "HospDeath", "scale")
    ) )   
}
 
# Ok to here

   needremoved=((as.character(data$observedType)=="med")&(as.character(data$type)=="M"))
   # depends if you're lost or not

if (sum(needremoved)>0)
   data[data$lost,"removed"]=round(rweibull(sum(data$lost),
  shape= getVecParams(params, "MedRec", "shape"),
    scale= getVecParams(params, "MedRec", "scale")
    ) )   
    

#if not lost
    needremoved=((as.character(data$observedType)=="med")&(as.character(data$type)=="M")&(data$lost==F))
if (sum(needremoved)>0)
     data[needremoved,"removed"]=censorweibull(
getVecParams(params, "MedRec", "shape"),
 getVecParams(params, "MedRec", "scale"),
(data[needremoved,"censor"]-data[needremoved,"med"])
    )   

   needremoved=((as.character(data$observedType)=="hosp")&(as.character(data$type)=="S"))

if (sum(needremoved)>0)
   data[needremoved,"removed"]=censorweibull(
getVecParams(params, "HospRec", "shape"),
 getVecParams(params, "HospRec", "scale"),
(data[needremoved,"censor"]-data[needremoved,"hospital"]-data[needremoved,"med"])
    ) 

   needremoved=((as.character(data$observedType)=="hosp")&(as.character(data$type)=="D"))

if (sum(needremoved)>0)
   data[needremoved,"removed"]=censorweibull(
getVecParams(params, "HospDeath", "shape"),
 getVecParams(params, "HospDeath", "scale"),
(data[needremoved,"censor"]-data[needremoved,"hospital"]-data[needremoved,"med"])
    ) 



data

}

