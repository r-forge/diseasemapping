rangeDistHazard=function(x, quantiles =c(0.025, 0.975), xseq=0:20) {

paramRange =list()
for(D in names(x)) {
 paramRange[[D]] = dprior(quantiles, x[[D]], "q")
  paramRange[[D]] =c(paramRange[[D]], x[[D]]["mean"])
  names(paramRange[[D]]) =c("lower", "upper", "mean")
}      


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

  # the prior mean
    theseParams = c(mean= paramRange$mean["mean"],
      shape= paramRange$shape["mean"], 
      scale = paramRange$mean["mean"] / 
        gamma(1 + 1/paramRange$shape["mean"])  )
     names(theseParams) = gsub("\\.[[:alpha:]]+$", "", names(theseParams))


distn[,"mean"] =  dweibullRound(xseq,  theseParams)
    hazard[,"mean"] =  distn[,"mean"] / 
      (1-pweibullRound(xseq,  theseParams))

    hazard[is.na(hazard)]=NA
     hazard[hazard==Inf]=NA

     return(list(dist=distn, hazard=hazard, paramRange=paramRange))
}       
       
plotPrior=function(x, posteriorSample = NULL, 
  file=NULL, quantiles =c(0.025, 0.975), tex=FALSE,
  transition=NULL) {
                                                           
if(any(names(x)==transition)) {
  x = x[[transition]]
}




includeGraphicsString = "\\includegraphics[width=0.3\\textwidth]{"
if(is.null(file)) {
    if(length(x)==3)
      par(mfrow=c(2,3))
    if(length(x)==4)
      par(mfrow=c(2,3))
}
     

if(tex) {
  cat("\\begin{figure} \n\n")

}

     
for(D in names(x)) {
Dvec = paste(transition, ".", D, sep="")
  thedist = attributes(x[[D]])$distribution
  if(is.null(posteriorSample)) {
  if(thedist == "gamma")
    xseq = seq(0, round(5*x[[D]]["mean"]), len=100)     
  if(thedist == "beta")
    xseq = seq(0, 1, len=100)      
  } else {
    therange = range(c(range(posteriorSample[,Dvec]), 0,x[[D]]["mean"]))
    xseq = seq(therange[1],  therange[2], len=100)
  }  
  
if(!is.null(file)) {
  filename = paste(file, D, ".pdf", sep="")  
  pdf(filename, width=5,height=3)
  par(mar=c(2,2,0,0))
}    
  if(is.null(posteriorSample)) {
  plot(xseq, dprior(xseq, x[[D]]), type="l", 
    ylab='prob', xlab=D, main="")
  } else {
    hist(posteriorSample[,Dvec], breaks=20,    
      ylab='prob', xlab=D, main="", xlim=therange,
      probability=T)

    lines(xseq, dprior(xseq, x[[D]]), col="red")
    
  }  
    
if(!is.null(file)) {
  dev.off()
}    

if(tex) {
  cat("\\subfigure[", D, 
    "]{", includeGraphicsString, filename, "} }\n" , sep="" )

}



} # end loop through names of x     



if(all(c("mean","shape")%in% names(x))) {
  xseq <- 0:20


if(is.null(posteriorSample)) {

  distHazard = rangeDistHazard(x, quantiles, xseq)
  hazard = distHazard$hazard
  distn= distHazard$dist
  paramRange= distHazard$paramRange
  
            
} else { # have posterior sample

if(is.null(transition)) {
  warning("not sure which transition to use from the posterior sample")  
}
 

# bivariate normal based posterior confidence region
themean = apply(posteriorSample[,c("HospDeath.mean","HospDeath.shape")],2,mean)
    thecontour = ellipse(
      var(posteriorSample[,c("HospDeath.mean","HospDeath.shape")]),
      centre=themean)
    Ncontour = dim(thecontour)[1]
    thecontour = thecontour[
      c(order(thecontour[,1])[c(1, Ncontour)],
        order(thecontour[,2])[c(1, Ncontour)]),]
    
    thecontour[thecontour<=0]=0.001
    
    hazard = distn = matrix(NA, length(xseq), 5, 
 dimnames = list(NULL, 
    c(paste("m=", thecontour[,1], ",s=", thecontour[,2], sep=""), "mean" ) )
    )
     
    
    for(Dcontour in seq(1, dim(thecontour)[1])) {
        thisCol = paste("m=",  thecontour[Dcontour,1], ",s=",
          thecontour[Dcontour,2], sep="")
          
      theseParams = c(mean= thecontour[Dcontour,1],
      shape= thecontour[Dcontour,2], 
      scale = thecontour[Dcontour,1] / 
        gamma(1 + 1/thecontour[Dcontour,2])  )

         names(theseParams) = gsub("\\.[[:alpha:]]+$", "", names(theseParams))
         names(theseParams) = gsub("\\.[[:alpha:]]+$", "", names(theseParams))

    distn[,thisCol] =  dweibullRound(xseq,  theseParams)
            
    hazard[,thisCol] =  distn[,thisCol] / 
      (1-pweibullRound(xseq,  theseParams))  
    }

    thisCol = "mean"
     theseParams = c(mean= themean[1],
      shape= themean[2], 
      scale = themean[1] / 
        gamma(1 + 1/themean[2])  )
         names(theseParams) = gsub("\\.[[:alpha:]]+$", "", names(theseParams))
         names(theseParams) = gsub("\\.[[:alpha:]]+$", "", names(theseParams))




      distn[,thisCol] =  dweibullRound(xseq,  theseParams)
            
    hazard[,thisCol] =  distn[,thisCol] / 
      (1-pweibullRound(xseq,  theseParams))  
  } # end if have posterior sample
} # end if mean and shape in x

# change any NaN or Inf's to NA
    hazard[is.na(hazard)]=NA
     hazard[hazard==Inf]=NA




if(!is.null(file)) {
filename<- paste(file, "dist.pdf", sep="")
  pdf(filename, width=5,height=4)
  par(mar=c(2,2,0,0))
}    



if(!all(is.na(distn[,"mean"]))){ 
    theylim = c(0, 1.2*max(distn[,"mean"], na.rm=T))
} else {
  theylim=c(0,max(distn, na.rm=T))
}

matplot(xseq, distn, lwd=c(1,1,1,1,2), type="l", 
  col=c("grey","yellow","green","orange","black"), lty=1,
  ylim=theylim)
if(!is.null(file)) {
  dev.off()
}    

if(tex) {
  cat("\\subfigure[distribution]{", includeGraphicsString,  filename, "} }\n", sep="" )

}

filename<- paste(file, "hazard.pdf", sep="")
if(!is.null(file)) {
  pdf(filename, width=5,height=4)
  par(mar=c(2,2,0,0))
}    

if(!all(is.na(hazard[,"mean"]))){ 
    theylim = c(0, 1.2*max(hazard[,"mean"], na.rm=T))
} else {
  theylim=c(0,quantile(hazard, 0.95, na.rm=T))
}


matplot(xseq, hazard, lwd=c(1,1,1,1,2), type="l" ,  
col=c("grey","yellow","green","orange","black"),lty=1,
ylim=theylim)
if(!is.null(file)) {
  dev.off()
}    
if(tex) {
  cat("\\subfigure[hazard]{", includeGraphicsString,  filename, "} }\n" , sep="" )

}


if(tex) {
  cat("\\caption{", file, "}\n",
  "\\label{fig:", file, "}\n",
    "\\end{figure}\n", sep="" )

}


return(invisible())


}
