                       # epidemic is growing by proposalOffset cases per day
updateBinomialPoisson <- function(priorMean, YobsToday, probObs, N, runs, nthin, proposalOffset) { # we're never updating!

   vecN <- NULL

   for (i in 1:runs) {

      newN <- N
      
     
      newN <- rpois(1, N + proposalOffset)  # proposal distribution

      # needs to be fixed, the code gets stopped here!
         while (newN < YobsToday) { 
            newN <- rpois(1, N)
            # print(newN)
         }
             
      like <- prod(dbinom(YobsToday, newN, probObs))/prod(dbinom(YobsToday, N, probObs))
      priorN <- dpois(newN, priorMean)/dpois(N, priorMean)  
      propose <- (dpois(newN, priorMean)/ppois(YobsToday - .1, priorMean, lower = F))/  # prob going from N to Nnew
         (dpois(N, newN +  proposalOffset)/ppois(YobsToday - .1, newN + proposalOffset, lower = F))    # prob going from Nnew to N

      pi <- priorN*propose # posterior distribution
      piPropose <- pi*propose
      alpha <- min(1, piPropose)

      if (alpha >= runif(1)) 
         {N <- newN}
      
      vecN <- c(vecN, N)
      
   }
   
   vecN[seq(from=nthin, to=runs,by=nthin)]

}
  
                                                                                                                                      
probSumWeibulls = function(parameters, x, Nsim) {

   params <- parameters[dim(parameters)[1],]  
   
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
	  
	  overX = rWeibullSum < x
	  
	  overSample[[Dtype]]] = cbind(InfOns=rWeibull1[x,], OnsMed=rWeibull2[x,],sum=rWeibullSum[x,])
	  
	  
      prob[Dtype] <- mean(overX)
   }

   list("prob" = prob, sample=overSample)
                                                             
}

probSumWeibulls(mcmcPandemicRun, 5, 1000)
    