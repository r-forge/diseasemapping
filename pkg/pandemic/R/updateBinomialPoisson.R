
                       # epidemic is growing by proposalOffset cases per day
updateBinomialPoisson <- function(priorMean, Y, prob, 
  runs, nthin=runs, proposalOffset=.2,Nstart=floor(Y / prob)) { # we're never updating!

  if (prob==1) {
   Yunique = unique(Y)
   if(length(Yunique)==1) {
    return(rep(Y, ceiling(runs/nthin)))
   } else {
      warning("prob=1 but different values of Y given")
   }
  }
  
  if (prob==0) {
  
      return(rpois(ceiling(runs/nthin), lambda = priorMean))     #return(rep(Inf,ceiling(runs/nthin)))
  
  }


   vecN <- NULL

   N = Nstart

   for (i in 1:runs) {

      newN <- N
      
      newN <- rpois(1, N + proposalOffset)  

         while (newN < Y-0.5) { 
            newN <- rpois(1, N + proposalOffset)
            # print(newN)
         }
             
      like <- exp(sum(dbinom(Y, newN, prob, log=T)) - sum(dbinom(Y, N, prob, log=T)) )       # problem because N > Y
      priorN <- dpois(newN, priorMean)/dpois(N, priorMean)  
      propose <- (dpois(newN, priorMean)/ppois(Y - .1, priorMean, lower = F))/  # prob going from N to Nnew
         (dpois(N, newN +  proposalOffset)/ppois(Y - .1, newN + proposalOffset, lower = F))    # prob going from Nnew to N

      PI <- priorN*like # posterior distribution
      PIpropose <- PI*propose
      alpha <- min(1, PIpropose)

    #  cat(newN, " ", like, " ", priorN, " ",  alpha, "\n")


      if (alpha >= runif(1)) 
         {N <- newN}
      
      vecN <- c(vecN, N)
      
   }
   
   vecN[seq(from=nthin, to=runs,by=nthin)]

}
  
                                                                                                                                      
probSumWeibulls = function(params, x, Nsim) {

 #  params <- parameters[dim(parameters)[1],]  
   
   prob = rep(NA, 3)
   names(prob) = c("M","S","D")

   pasteShapeTransition1 <- paste("InfOns", "shape", sep = ".")
   pasteScaleTransition1 <- paste("InfOns", "scale", sep = ".")

   rWeibull1 <- rweibull(Nsim, shape = params[pasteShapeTransition1], scale = params[pasteScaleTransition1])
   
   overSample = list()
   
   for(Dtype in names(prob)) {

      trans = paste("OnsMed", Dtype, sep="")
   
      pasteShapeTransition2 <- paste(trans, "shape", sep = ".")
      pasteScaleTransition2 <- paste(trans, "scale", sep = ".")

      rWeibull2 <- rweibull(Nsim, shape = params[pasteShapeTransition2], scale = params[pasteScaleTransition2])

	rWeibullSum = rWeibull1 + rWeibull2
	  
	  overX = rWeibullSum > x
	  
	  overSample[[Dtype]] = cbind(InfOns=rWeibull1[overX], OnsMed=rWeibull2[overX], sum=rWeibullSum[overX])
	 
      prob[Dtype] <- mean(overX)
	rWeibullSum = rWeibull1 + rWeibull2
	  
	  overX = rWeibullSum < x
	  
	  overSample[[Dtype]]] = cbind(InfOns=rWeibull1[x,], OnsMed=rWeibull2[x,],sum=rWeibullSum[x,])
	  
	  
      prob[Dtype] <- mean(overX)

   }


   list("prob" = prob, "sample" = overSample)

                                                            
}

              
              
              
              
              
              
              
              
              
              
              
              
              
              
              
              
              
              
              
              
              
              
              
              
