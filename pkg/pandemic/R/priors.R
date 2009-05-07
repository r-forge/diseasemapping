                                        pandemicPriors = function( 
  InfOns = meanShapeZerosPrior(),
 OnsMedM = meanShapeZerosPrior(),
 OnsMedS = meanShapeZerosPrior(),
 OnsMedD = meanShapeZerosPrior(),
 MedRec = meanShapeZerosLostPrior(),
 MedHospS = meanShapeZerosPrior(),
 MedHospD = meanShapeZerosPrior(),
 HospRec = meanShapeZerosPrior(),
 HospDeath = meanShapeZerosPrior(),
 probs = probsPrior()
)
{
  theFormals = formals()
  result = list()
  for(D in names(theFormals))
     result[[D]] = eval(theFormals[[D]])
  result   
}


meanShapeZerosPrior = function(
  mean =  c(mean=1, sd=1),
  shape = c(mean=1, sd=1),
  zeros =  c(mean=1, sd=1)
) {

  theFormals = formals()
  result = list()
  for(D in names(theFormals)) {
     result[[D]] = eval(theFormals[[D]])
   }
   priorShapeScale(result)
}

 meanShapeZerosLostPrior = function(lost=c(mean=1, sd=1), ...) {
  list(???)
 
 }


priorShapeScale = function(priorList) {
# prior for one transition, priorList has elements mean, shape, zeros, lost
# each element has "mean" and "sd"

  thenames = names(priorList)
  for(D in thenames[thenames %in% c("mean","shape") ) {
    priorList[[D]]["shape"] = priorList[[D]]["sd"]^2/priorList[[D]]["mean"]
    priorList[[D]]["scale"] = priorList[[D]]["mean"] /priorList[[D]]["shape"]
  }
  for(D in thenames[thenames %in% c("zeros","lost") ) {
     priorList[[D]]["shape2"] = 
      priorList[[D]]["sd"]^2 / priorList[[D]]["mean"] / (1-priorList[[D]]["mean"])*
        (1+1/(1-priorList[[D]]["mean"]))
     priorList[[D]]["shape1"] = priorList[[D]]["shape2"] * 
      priorList[[D]]["mean"] / (1-priorList[[D]]["mean"])
  }
  priorList
}


