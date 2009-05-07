pandemicParams <- function(
 InfOns = c(mean=1,shape = 1, zeros = 0),
 OnsMedM = c(mean = 1, shape = 1, zeros = 0),
 OnsMedS = c(mean = 1, shape = 1, zeros = 0),
 OnsMedD = c(mean = 1, shape = 1, zeros = 0),
 MedRec = c(mean = 1, shape = 1, zeros = 0),
 MedHospS = c(mean = 1, shape = 1, zeros = 0),
 MedHospD = c(mean = 1, shape = 1, zeros = 0),
 HospRec = c(mean = 1, shape = 1, zeros = 0),
 HospDeath = c(mean = 1, shape = 1, zeros = 0)
) {
  theFormals = formals()
  result = list()
  for(D in names(theFormals))
     result[[D]] = eval(theFormals[[D]])
   result               

}


addScaleParameters = function(params) {

  for(D in names(params))
    params[[D]]["scale"] = params[[D]]["mean"] / 
      gamma(1 + 1/params[[D]]["shape"])
  params
}

addMeanParameters = function(params) {

  for(D in names(params))
    params[[D]]["mean"] = params[[D]]["scale"] *
      gamma(1 + 1/params[[D]]["shape"])
  params
}