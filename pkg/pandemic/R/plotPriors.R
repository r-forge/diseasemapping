       
plotPrior=function(x,file=NULL, quantiles =c(0.025, 0.975)) {
if(is.null(file)) {
    if(length(x)==3)
      par(mfrow=c(2,3))
    if(length(x)==4)
      par(mfrow=c(2,2))
} else {
  
}     

params = list()
     
for(D in names(x)) {
  thedist = attributes(x[[D]])$distribution
  if(thedist == "gamma")
    xseq = seq(0, 5*x[[D]]["mean"], len=100)     
  if(thedist == "beta")
    xseq = seq(0, 1, len=100)      
  plot(xseq, dprior(xseq, x[[D]]), type="l", main=D,
    ylab='prob', xlab='value')

 # get conf int for parameters
paramRange[[D]] = dprior(quantiles, x[[D]], "q")
paramRange[[D]] =c(paramRange[[D]], x[[D]]["mean"])
names(paramRange[[D]]) =c("lower", "upper", "mean")

}      



if(all(c("mean","shape")%in% names(paramRange))) {
hazard = distn = matrix(NA, length(xseq), 5, 
 dimnames = list(NULL, 
    c(outer(paste("m=",paramRange$mean[c("lower", "upper")] ,sep=""), 
        paste("s=",paramRange$shape[c("lower", "upper")] ,sep=""),
        FUN=paste), "mean") ) )

for(Dmean in c("upper","lower")) {


}

return(distn)

}




# hazard function

}