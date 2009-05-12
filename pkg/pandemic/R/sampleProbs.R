
sampleProbs <- function(latentTypes, probPrior) {

# fatality
latentTypes = latentTypes[!is.na(latentTypes)]

fatal = latentTypes=="D"
Nfatal = sum(fatal)
fatal = c("TRUE"=Nfatal, "FALSE"=length(fatal) - Nfatal)

fatal= rbeta(1,
shape1=probPrior$fatality["shape1"] + fatal["TRUE"],
shape2=probPrior$fatality["shape2"] + fatal["FALSE"])



# hospital
hosp = latentTypes[latentTypes!= "D"] == "S"
Nhosp = sum(hosp)
hosp = c("TRUE"=Nhosp, "FALSE"=length(hosp)-Nhosp)


hospGivenAlive = rbeta(1,
shape1=probPrior$hosp["shape1"] + hosp["TRUE"],
shape2=probPrior$hosp["shape2"] + hosp["FALSE"]
)

c(M=(1-fatal)*(1-hospGivenAlive), S=(1-fatal)*hospGivenAlive,
  D=fatal)

}