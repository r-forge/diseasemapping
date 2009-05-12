pandemicParams <- function(
 InfOns = c(mean=1,shape = 1, zeros = 0),
 OnsMedM = c(mean = 3, shape = 2, zeros = 0),
 OnsMedS = c(mean = 2.5, shape = 1, zeros = 0),
 OnsMedD = c(mean = 2, shape = 1, zeros = 0),
 MedRec = c(mean = 7, shape = 2, zeros = 0, lost=0.9),
 MedHospS = c(mean = 2.5, shape = 2, zeros = 0),
 MedHospD = c(mean = 1.5, shape = 1.5, zeros = 0),
 HospRec = c(mean = 4, shape = 1.5, zeros = 0),
 HospDeath = c(mean = 3, shape = 1, zeros = 0),
 probs = c(M=0.6, S=0.3, D=0.1)
) {
  theFormals = ls(-1)
  result = list()
  for(D in theFormals)
     result[[D]] = get(D, pos=-1)
   
   addScaleParameters(result)               

}




addScaleParameters = function(params) {

  thenames = names(params)
  thenames = thenames[thenames != "probs"]
  for(D in thenames)
    params[[D]]["scale"] = params[[D]]["mean"] / 
      gamma(1 + 1/params[[D]]["shape"])
  params
}

addMeanParameters = function(params) {

  thenames = names(params)
  thenames = thenames[thenames != "probs"]
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