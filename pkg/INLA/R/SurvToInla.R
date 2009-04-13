SurvToInla = function(x) {

if(class(x) != "Surv")
  warning("x is of class ", class(x), ", should be Surv")

Ncol = dim(x)[2]
Nobs = dim(x)[1]
type = attributes(x)$type

NAvec = rep(-1, Nobs)
result <- data.frame(truncation=rep(0,Nobs),
  lower=NAvec, upper=NAvec, event=NAvec, 
  time = NAvec) 

xevent <- x[,"status"]
hadEvent = xevent==1
result$event = xevent
  
if(type != "counting")
  result[hadEvent, "time"] = x[hadEvent,"time"]

if(Ncol==2) {
 # either right or left censored

   if(type=="left") {
     result[!hadEvent,"event"] = 2
     result[!hadEvent,"upper"]= x[!hadEvent,"time"]
   } else if(type=="right") {
     result[!hadEvent,"lower"]= x[!hadEvent,"time"]
   } else {
      warning("type is ", type, ", should be left or right")
   }
} else if(Ncol==3) {
# either interval censored or left truncated
   result$event = xevent
   if(type=="interval") {
      # rights
        rights = xevent==0
        result[rights,"lower"]=x[rights,"time"]
      # lefts
        lefts  = xevent==2
        result[lefts,"upper"]=x[lefts,"time"]
      # intervals
        intervals = xevent==3
        result[intervals,"lower"]=x[intervals,"time"]
        result[intervals,"upper"]=x[intervals,2]
         
   } else if(type=="counting") {
        result$truncation = x[,1]
        result$event = xevent
        result[hadEvent, "time"] = x[hadEvent,2]
        result[!hadEvent, "lower"] = x[hadEvent,2]
        
   } else {
      warning("type is ", type, ", should be interval or counting")
   }

} else {
  warning(Ncol," columns in x, should be 2 or 3")
}

return(result)
}
