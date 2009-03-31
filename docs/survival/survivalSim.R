  
	n = 10000
	nu = 1.5

	
simData = function(n, nu, nscale=2, truncate=T)	{
# if truncate=T, left truncated right censored data are produced
# if F, then interval censored data result.
  
  if(truncate)
    nbig = n*nscale 
  else
    nbig = n
    
	z = runif(nbig)
	gamma = exp(1 + z)
	theta = gamma^(-nu)


  result = data.frame(time = rweibull(nbig, nu, theta),
    z=z)
    
  if(truncate) {    
  threashold = quantile(result$time, 0.5)
  
  result$Li = runif(nbig, 0, threashold)
  followup = runif(nbig, threashold, 
    quantile(result$time, 0.95))
   
  result$Ri = pmin(followup, result$time)

  result$Event = result$time < followup  

  # remove observations which have events before Li
  result = result[result$Li < result$time,]

  result$Event = as.integer(result$Event)
  result$typeOfEvent = factor(result$Event,
    levels=c(1,0), labels=c("left truncated and observed event", "left truncated and right censored"))  
  } else {
    RightCensored = sample(nbig, round(nbig*0.6))
    IntervalCensored = sample(RightCensored, 
      round(length(RightCensored)*.4))
    result$Event = 1
    result$Event[RightCensored]=0
    result$Event[IntervalCensored]=3
    result$typeOfEvent = factor(result$Event,
      levels=c(0,1,3), labels=c("right censored", "observed event", "interval censored"))  
    
    result$Li = NA
    result$Ri = result$time
    result$Li[IntervalCensored] = runif(length(IntervalCensored),
      0, result$time[IntervalCensored])
    result$Ri[RightCensored] = runif(length(RightCensored),
      result$time[RightCensored], max(result$time))
  }
  
  if(n > dim(result)[1])
    print("not enough non-truncated observations, make nscale bigger")
  
  result = result[1:n,]

      

result
}

x = simData(n, nu)
write.table(x, "survTruncated.txt")
xInt = simData(n, nu, truncate=F)
write.table(xInt, "survInterval.txt")
