pandemicParams <- function(
 InfOns = c(mean=1,shape = 1, zeros = 0.1),
 OnsMedM = c(mean = 3, shape = 2, zeros = 0.1),
 OnsMedS = c(mean = 2.5, shape = 1, zeros = 0.1),
 OnsMedD = c(mean = 2, shape = 1, zeros = 0.1),
 MedRec = c(mean = 7, shape = 2, zeros = 0.1, lost=0.9),
 MedHospS = c(mean = 2.5, shape = 2, zeros = 0.1),
 MedHospD = c(mean = 1.5, shape = 1.5, zeros = 0.1),
 HospRec = c(mean = 4, shape = 1.5, zeros = 0.1),
 HospDeath = c(mean = 3, shape = 1, zeros = 0.1),
 probs = c(M=0.6, S=0.3, D=0.1),
 ageProbs=NULL
) {
  theFormals = ls(-1)
  result = list()
  for(D in theFormals)
     result[[D]] = get(D, pos=-1)
  
  if(!all(names(ageProbs) %in% c("S","D")) )
  {
    cat(names(ageProbs))
   warning(" names of ageProbs must be either S or D")
  }
   
   addScaleParameters(result)               

}

addAgeProbs =function(age=0:100, prob=rep(0.1, length(age)) ) {
     data.frame(age, prob)
}


addScaleParameters = function(params) {

  thenames = names(params)
  thenames = thenames[-grep("(age)?[pP]robs$", thenames)]
  for(D in thenames)
    params[[D]]["scale"] = params[[D]]["mean"] / 
      gamma(1 + 1/params[[D]]["shape"])
  params
}

addMeanParameters = function(params) {

  thenames = c("probs",names(params))
  thenames = thenames[-grep("(age)?[pP]robs$", thenames)]
  for(D in thenames)
    params[[D]]["mean"] = params[[D]]["scale"] *
      gamma(1 + 1/params[[D]]["shape"])
  params
}

getVecParams = function(params, string="OnsMed", value="shape") {
  result = unlist(
    lapply(params[grep(paste("^", string, sep=""), names(params))], 
      function(qq) qq[value])
    )
  names(result) = gsub(paste("^", string, "|.", value, "$|.NA", sep=""), 
    "", names(result))
  result  


}


vecParamsToList = function(vecParams) {
  thenames = names(vecParams)
  
  ageProbsIndex = grep("^ageProbs", thenames)
  if(length(ageProbsIndex)) {
    namesNotAgeProbs = thenames[-ageProbsIndex]
    namesAgeProbs = thenames[ageProbsIndex]
  } else {
   namesNotAgeProbs = thenames
   namesAgeProbs=NULL
  }
  
  matNames = matrix(unlist(strsplit(namesNotAgeProbs, "\\.")) , ncol=2, byrow=T)
  
  Stransition = unique(matNames[,1])
  
  params = list()
  for(D in Stransition) {
  Dstring <- paste("^", D, "\\.", sep="")
    forList = vecParams[grep(Dstring, namesNotAgeProbs, value=T)]
    names(forList) = gsub(Dstring, "", names(forList))
    params[[D]] = forList
  }
  
  if(length(namesAgeProbs)) {
    for(Dtrans in c("S", "D") ){
      thisTrans = grep(paste("^ageProbs.",Dtrans,".prob[[:digit:]]+", sep=""),
        namesAgeProbs, value=T)
      thisAge = as.numeric(gsub( paste("^ageProbs.",Dtrans,".prob", sep=""), 
        "", thisTrans))
      theorder = order(thisAge)
      params$ageProbs[[Dtrans]] = data.frame(
       age = thisAge[theorder],
       prob = vecParams[thisTrans[theorder]]
      )  
    }
  }
  params
}
