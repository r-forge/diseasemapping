simHospitals = function(paramSim, Ndays=100, NinfectionsPerDay=10) {

  Nsim = dim(paramSim)[1]
  result = matrix(NA, Ndays,Nsim)
  for(Dsim in 1:Nsim) {
    params = vecParamsToList(paramSim[Dsim,])
    data = simEpidemic(params, delta=NinfectionsPerDay, days=Ndays,
      probOnsetMissing=0, randomInfections=T)
    data = data[data$observedType %in% c("D","S","hosp"),] 
    data$removed[is.na(data$removed)] = Inf
      
    data = data.frame(start=data$hospital+data$med, end=data$removed+data$med)
    data$end[is.na(data$end)] = Inf
    for(Dday in 1:Ndays) {
      notOut = data$end >= Dday
        isIn = data$start <= Dday
      result[Dday,Dsim] = sum(notOut & isIn)
      data = data[notOut,]
    }
  
  }  

  result = list(sample = result,
    summary = t(apply(result, 1, function(qq) {
      c(mean=mean(qq), quantile(qq, probs=c(0.025, 0.5, 0.975)))
    } ) )
    )
}


plotHospitals = function(x) {
  if(is.list(x)) x = x$summary
  days = seq(1, dim(x)[1])
  
  plot(NA, xlim=range(days), ylim=range(x), xlab="days", ylab="beds")
  
  polygon(c(days, rev(days)), c(x[,"2.5%"], rev(x[,"97.5%"])),
    col="yellow", lty=0)
  lines(days, x[,"mean"])  
  
}