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
    # start = day (of study) individual went into hospital
    # end = day (of study) individual went into the hospital
    data$end[is.na(data$end)] = Inf
    for(Dday in 1:Ndays) {
      notOut = data$end >= Dday
      isIn = data$start <= Dday
    # if notOut is TRUE, individual is still in the hopital; if FALSE, individual is out of the hospital
    # if isIN is TRUE, individual is in the hospital; if FALSE, individual is not yet in the hospital
    # this is recorded for every one of Ndays (the number of research days) 
      result[Dday,Dsim] = sum(notOut & isIn)   
    # in the result matrix, the (Dday, Dsim) entry gets a value of 1 if notOut and isIn are both TRUE and 0, otherwise 
    # (which means that individual has entered the hospital and has not left yet)
    # but, we sum along the columns, so where the column is 1, individual enters, the final number is the total number of individuals in the hospital on that day of the study
      data = data[notOut,]   # output data, where individual is not out yet out of the hospital
    }
  
  }  

  result = list(sample = result,
    summary = t(apply(result, 1, function(qq) {
      c(mean=mean(qq), quantile(qq, probs=c(0.025, 0.5, 0.975)))  # return the mean, 0.025, 0.5, 0.975 quantiles along the rows of the result matrix
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