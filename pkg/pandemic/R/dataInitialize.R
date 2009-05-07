censorweibull=function(a,b,cen)
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

   data[needremoved,"removed"]=censorweibull(
getVecParams(params, "MedRec", "shape"),
 getVecParams(params, "MedRec", "scale"),(data[needremoved,"censor"]-data[needremoved,"med"])
    ) 


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




