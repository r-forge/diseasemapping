plotAgeProbs <- function(postSample, type="fatality", 
  quantiles=c(0.025, 0.975),prior=NULL) {

  typeSD = c(fatality="D", hosp="S")[type]
  
  postSample = postSample[,
    grep(paste("^ageProbs\\.", typeSD, "\\.prob[[:digit:]]+$", sep=""),
      colnames(postSample))] 
  
  age = as.numeric(
    gsub(paste("^ageProbs\\.", typeSD, "\\.prob",sep=""), "", colnames(postSample) )
   )
   
  toplot = as.data.frame(t(apply(postSample, 2, quantile, probs=quantiles) ))
  names(toplot) = gsub("\\%$", "pct", names(toplot))
  
  toplot$mean = apply(postSample, 2, mean)
  
  matplot(age, toplot, type="n", xlab="age", 
    ylab =paste("prob(", type, ")")) 
  matlines(age,  t(postSample), col="grey", lwd=0.5)
  matlines(age, toplot,
    col="black", lty=c(2,2,1)) 
    
}
