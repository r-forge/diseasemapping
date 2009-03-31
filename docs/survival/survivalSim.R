  if(F) {
  
	n = 10000
	nu = 1.5
 x = simData(n, nu)
write.table(x, "survTruncated.txt")
xInt = simData(n, nu, truncate=F)
write.table(xInt, "survInterval.txt")

library(survival)
xInt = simData(30000, 1.5, truncate=F)
# Surv wants right censoring times in time, not time2
#xInt$Li[is.na(xInt$Li)] = xInt$Ri[is.na(xInt$Li)]
survdata = Surv(xInt$L, xInt$R, xInt$Event, type="interval")
xnoI = xInt[xInt$Event !=0,]
survnoI = Surv(xnoI$Li, xnoI$R, xnoI$Event, type="interval")

# offests for survreg, because of different parametrization
# multiply by -1/nu
zSurvReg = (-1/nu) * xInt$z 
zSurvRegnoI = (-1/nu) * xnoI$z 
	
res=survreg(formula = survdata ~ 1 + offset(zSurvReg))
resnoI=survreg(formula = survnoI ~ 1 + offset(zSurvRegnoI))

# this should be 1
-res$coef/res$scale
-resnoI$coef/resnoI$scale

# this should be 1.5
1/res$scale
1/resnoI$scale

}	
	
simData = function(n, nu, nscale=2, truncate=T)	{
# if truncate=T, left truncated right censored data are produced
# if F, then interval censored data result.
  
  if(truncate)
    nbig = n*nscale 
  else
    nbig = n
    
	z = runif(nbig)
	gamma = exp(1 + z)
	theta = gamma^(-1/nu)


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
   
  if(n > dim(result)[1])
    warning("not enough non-truncated observations, make nscale bigger")
  
  result = result[1:n,]

  result$Event = as.integer(result$Event)
  result$typeOfEvent = factor(result$Event,
    levels=c(1,0), labels=c("left truncated and observed event", "left truncated and right censored"))  
  } else {   # interval censored
    result$Event = 1

    RightCensored = sample(nbig, round(nbig*0.6))
    IntervalCensored = sample(RightCensored, 
      round(length(RightCensored)*.7))
    result$Event = 1
    result$Event[RightCensored]=0
    result$Event[IntervalCensored]=3
    result$typeOfEvent = factor(result$Event,
      levels=c(0,1,3), labels=c("right censored", "observed event", "interval censored"))  
    
    result$Li = result$time
    result$Li[RightCensored] = pmax(
      floor(result$time[RightCensored]*10)/10,
      result$time[RightCensored]/100)

    result$Ri = NA
    result$Ri[IntervalCensored] = 
      ceiling(result$time[IntervalCensored]*10)/10

  }


      

result
}


