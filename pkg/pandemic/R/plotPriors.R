       
plotPrior=function(x,file=NULL, quantiles =c(0.025, 0.975), tex=FALSE) {
if(is.null(file)) {
    if(length(x)==3)
      par(mfrow=c(2,3))
    if(length(x)==4)
      par(mfrow=c(2,2))
} else {
  
}     

if(tex) {
  cat("\\begin{figure} \n\n")

}

params = paramRange = list()
     
for(D in names(x)) {
  thedist = attributes(x[[D]])$distribution
  if(thedist == "gamma")
    xseq = seq(0, round(10*x[[D]]["mean"]), len=100)     
  if(thedist == "beta")
    xseq = seq(0, 1, len=100)      
  
filename = paste(file, D, ".pdf", sep="")  
if(!is.null(file)) {
  pdf(filename, width=5,height=3)
  par(mar=c(2,2,0,0))
}    
  plot(xseq, dprior(xseq, x[[D]]), type="l", 
    ylab='prob', xlab=D, main="")
if(!is.null(file)) {
  dev.off()
}    

if(tex) {
  cat("\\subfigure[", D, 
    "]{ \\includegraphics{", filename, "} }\n" )

}



 # get conf int for parameters
paramRange[[D]] = dprior(quantiles, x[[D]], "q")
paramRange[[D]] =c(paramRange[[D]], x[[D]]["mean"])
names(paramRange[[D]]) =c("lower", "upper", "mean")

}      




if(all(c("mean","shape")%in% names(paramRange))) {
xseq <- 0:20

hazard = distn = matrix(NA, length(xseq), 5, 
 dimnames = list(NULL, 
    c(outer(paste("m=",paramRange$mean[c("lower", "upper")] ,sep=""), 
        paste("s=",paramRange$shape[c("lower", "upper")] ,sep=""),
        FUN=paste), "mean") ) )

for(Dmean in c("upper","lower")) {
  for(Dshape in c("upper","lower")) {
    thisCol = paste("m=",  paramRange$mean[Dmean], " s=",
          paramRange$shape[Dshape], sep="")
    theseParams = c(mean= paramRange$mean[Dmean],
      shape= paramRange$shape[Dshape], 
      scale = paramRange$mean[Dmean] / 
        gamma(1 + 1/paramRange$shape[Dshape])  )
      # get rid of the names suffixes which R adds for some reason  
     names(theseParams) = gsub("\\.[[:alpha:]]+$", "", names(theseParams))
    distn[,thisCol] =  dweibullRound(xseq,  theseParams)
            
    hazard[,thisCol] =  distn[,thisCol] / 
      (1-pweibullRound(xseq,  theseParams))


  }
}

 theseParams = c(mean= paramRange$mean["mean"],
      shape= paramRange$shape["mean"], 
      scale = paramRange$mean["mean"] / 
        gamma(1 + 1/paramRange$shape["mean"])  )
names(theseParams) = gsub("\\.[[:alpha:]]+$", "", names(theseParams))

distn[,"mean"] = dweibullRound(xseq,  theseParams)
            
hazard[,"mean"] =  distn[,"mean"] / 
      (1-pweibullRound(xseq,  theseParams))

hazard[hazard==Inf] = NA

filename<- paste(file, "dist.pdf", sep="")
if(!is.null(file)) {
  pdf(filename)
  par(mar=c(2,2,0,0))
}    

matplot(xseq, distn, lwd=c(1,1,1,1,2), type="l", 
  col=c("grey","yellow","green","orange","black"), lty=1,
  ylim=c(0, max(distn[,"mean"])))
if(!is.null(file)) {
  dev.off()
}    

if(tex) {
  cat("\\subfigure[distribution ]{ \\includegraphics{", filename, "} }\n" )

}

filename<- paste(file, "hazard.pdf", sep="")
if(!is.null(file)) {
  pdf(filename)
  par(mar=c(2,2,0,0))
}    

matplot(xseq, hazard, lwd=c(1,1,1,1,2), type="l" , col="black", lty=1,
ylim=c(0, max(hazard[,"mean"])))
if(!is.null(file)) {
  dev.off()
}    
if(tex) {
  cat("\\subfigure[hazard]{ \\includegraphics{", filename, "} }\n" )

}


#return(list(x=xseq, dist=distn, haz=hazard, par=theseParams))

}

if(tex) {
  cat("\\caption{", file, "}\n",
  "\\label{fig:", file, "}\n",
    "\\end{figure}\n", sep="" )

}


# hazard function

}