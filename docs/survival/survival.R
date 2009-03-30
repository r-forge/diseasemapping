  
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
  result$followup = runif(nbig, threashold, 
    quantile(result$time, 0.95))
  
  
  # remove observations which have events before Li
  result = result[result$Li < result$time,]
  result$Ri = pmin(result$followup, result$time)
  result$Event = result$time < result$followup  
  
  } else {
    censored = as.logical(rbinom(n, 1,0.5))
    result$Event = c(0,3)[1+censored]
    result$Li = NA
    result$Ri = result$time
    result$Li[censored] = runif(sum(censored),
      0, result$time)
    result$Ri[censored] = runif(sum(censored),
      result$time, max(result$time))
  }
  
  if(n > dim(result)[1])
    print("not enough non-truncated observations, make nscale bigger")
  
  result = result[1:n,]

result
}

x = simData(n, nu)
write.table(x, "surv.txt")
xInt = simData(n, nu, truncate=F)
write.table(xInt, "survInterval.txt")
