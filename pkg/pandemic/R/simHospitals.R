simHospitals = function(paramSim, Ndays=100, NinfectionsPerDay=100) {

  result = NULL
  Nsim = dim(paramSim)[1]
  for(D in 1:Nsim) {
    params = vecParamsToList(paramSim[Dsim,])
  data = simEpidemic(params, N, days=Ndays,
    probOnsetMissing=0)  
  
  
  }  


}