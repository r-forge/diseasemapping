censorweibull=function(a,b,cen, params)
{
weibull=0
for(i in 1:length(cen))
{
weibull[i]=0
while(weibull[i]<=cen[i])
{
weibull[i]=round(rweibull(1,a,b))
}
}
weibull
}




# hosps
inHosp = which(data$observedType == "hosp")
inHospTimes = data[inHosp, "censor"] - data[inHosp, "hospital"]
probCensorGivenDeadly = 1-pweibullRound(inHospTimes, params$HospDeath)
probCensorGivenSerious = 1-pweibullRound(inHospTimes, params$HospRec)

probDeadlyGivenCensor = probCensorGivenDeadly *params$probs["D"]

probDeadlyGivenCensor = probDeadlyGivenCensor / 
  (probDeadlyGivenCensor + probCensorGivenSerious)
  
data[inHosp,"type"] = c("S","D")[1+rbinom(length(inHosp), 1, probDeadlyGivenCensor)]  


# meds
inMed = which(data$observedType == "med")
inMedTimes = data[inMed, "censor"] - data[inMed, "med"]

probOf = data.frame(
  "L" = 1* params$probs["M"] * params$MedRec["lost"], 
  "M"= (1-pweibullRound(inMedTimes, params$MedRec)) *
      (params$probs["M"] * (1- params$MedRec["lost"]) ) ,
  "S"= params$probs["S"]*
    (1-pweibullRound(inMedTimes, params$MedHospS)),
  "D"= params$probs["D"]*
    (1-pweibullRound(inMedTimes, params$MedHospD) )
  )
#sumProb = apply(probOf, 1, sum)  

#probOf = apply(probOf, 2, function(qq) qq/sumProb)
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

   data[needhospital,"hospital"]=censorweibull(
    getVecParams(params, "MedHospS", "shape"),
    getVecParams(params, "MedHospS", "scale"),
(data[needhospital,"censor"]-data[needhospital,"med"])
    ) 

   
  data[needhospital,"removed"]  =  round(rweibull(sum(needhospital),
  shape= getVecParams(params, "HospRec", "shape"),
    scale= getVecParams(params, "HospRec", "scale")
    ) )   

   needhospital=((as.character(data$observedType)=="med")&(as.character(data$type)=="D"))

   data[needhospital,"hospital"]=censorweibull(
     getVecParams(params, "MedHospD", "shape"),
    getVecParams(params, "MedHospD", "scale"),
(data[needhospital,"censor"]-data[needhospital,"med"])
    ) 

   
  data[needhospital,"removed"]  =  round(rweibull(sum(needhospital),
  shape= getVecParams(params, "HospDeath", "shape"),
    scale= getVecParams(params, "HospDeath", "scale")
    ) )   


   needremoved=((as.character(data$observedType)=="med")&(as.character(data$type)=="M"))
   # depends if you're lost or not

   # if lost
   data[data$lost,"removed"]=round(rweibull(sum(data$lost),
  shape= getVecParams(params, "MedRec", "shape"),
    scale= getVecParams(params, "MedRec", "scale")
    ) )   
    #if not lost



   needremoved=((as.character(data$observedType)=="hosp")&(as.character(data$type)=="S"))

   data[needremoved,"removed"]=censorweibull(
getVecParams(params, "HospRec", "shape"),
 getVecParams(params, "HospRec", "scale"),
(data[needremoved,"censor"]-data[needremoved,"hospital"]-data[needremoved,"med"])
    ) 

   needremoved=((as.character(data$observedType)=="hosp")&(as.character(data$type)=="D"))

   data[needremoved,"removed"]=censorweibull(
getVecParams(params, "HospDeath", "shape"),
 getVecParams(params, "HospDeath", "scale"),
(data[needremoved,"censor"]-data[needremoved,"hospital"]-data[needremoved,"med"])
    ) 

data




