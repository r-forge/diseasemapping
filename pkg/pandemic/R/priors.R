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
    theFormals = ls(-1)
  result = list()
  for(D in theFormals)
     result[[D]] = get(D, pos=-1)
# check that the probs have a suitable prior
if(!all(c("fatality","hosp") %in% names(result$probs)))
  warning("prior for probabilities needs to have 'fatality' and 'hosp' elements")


  result   


}


meanShapeZerosPrior = function(
  mean =  c(mean=1, sd=1),
  shape = c(mean=1, sd=1),
  zeros =  c(mean=0.1, sd=0.1)
) {

  theFormals = ls(-1)
  result = list()
  for(D in theFormals) {
     result[[D]] = get(D, pos=-1)
   }
   priorShapeScale(result)
}

 meanShapeZerosLostPrior = function(lost=c(mean=0.5, sd=0.1), ...) {
  
  result = meanShapeZerosPrior(...)
  result$lost = lost
  priorShapeScale(result)
  
 }


priorShapeScale = function(priorList) {
# prior for one transition, priorList has elements mean, shape, zeros, lost
# each element has "mean" and "sd"

  thenames = names(priorList)
  for(D in thenames[thenames %in% c("mean","shape")] ) {
    priorList[[D]]["scale"] = priorList[[D]]["sd"]^2/priorList[[D]]["mean"]
    priorList[[D]]["shape"] = priorList[[D]]["mean"] /priorList[[D]]["scale"]
    attributes(priorList[[D]])$distribution = "gamma"
  }
  for(D in 
    thenames[thenames %in% c("zeros","lost","fatality","hosp")] ) {
    Var = priorList[[D]]["sd"]^2 
    mu = priorList[[D]]["mean"]
    mu1 = 1-mu
    
    
    priorList[[D]]["shape2"] = (mu1/Var) * (mu * mu1 - Var)
     priorList[[D]]["shape1"] = priorList[[D]]["shape2"] * 
      (mu/mu1)
    attributes(priorList[[D]])$distribution = "beta"      
  }
  
  
  priorList
}


probsPrior = function(
  fatality = c(mean=0.1, sd=0.2),
  hosp = c(mean=0.2, sd=0.2)
  ) {
  
  result = list(
    fatality=fatality,
    hosp=hosp
  )

  priorShapeScale(result)  
}

psProbPriors = function(
  fatality=psPrior(),
  hosp=psPrior()
){
  result = list(fatality=fatality, hosp=hosp)
}

psPrior = function(
 taub1=2,
 taub2=0.05,
 priorMean=0.01,
 upper95 = 0.02) {

 # COMMENTS FROM PSGAM DOCUMENTATION BELOW
     # mean =taub1 / taub2, var propto taub1/taub2^2
 
#               taub1, taub2:  reals giving the hyperparameters of the prior 
#                       distribution for the inverse of the
#                       variances, 1/sigmab ~ Gamma(taub1/2,taub2/2).
 #rgamma(alpha,beta)
#c=======================================================================                  
#c     This function generates a random gamma value.
#c     The parametrization is such that E(X)=alpha/beta  
# I think taub1=alpha is the shape, taub2=beta is range = 1/scale


 result = list(taub1=as.numeric(taub1), taub2=as.numeric(taub2), 
 priorMean=as.numeric(priorMean), upper95=as.numeric(upper95), 
  beta0=as.numeric(log(priorMean / (1-priorMean)) ))
 
 upperLogit = log(upper95/(1-upper95))

 result$Sbeta0 = ( (upperLogit - result$beta0)/2  ) 
  
   attributes(result)$distribution = "psPrior"
   result
}   


dprior = function(x, prior, prefix="d", ...) {
  argList = list(
    gamma=c("shape","scale"),
    beta = c("shape1", "shape2")
  )[[attributes(prior)$distribution]]
  
  theFun = get(paste(prefix, attributes(prior)$distribution, sep=""))
  
  toCall = list(theFun, x, ...)
  for(D in argList)
    toCall[[D]] = prior[D]

  eval(as.call(toCall  ))
       
    
  
}