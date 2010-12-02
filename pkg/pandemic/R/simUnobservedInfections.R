

simUnobservedInfections = function(postSample, lengthOfEpidemic, infectParams = c(rate = 1, immigration = 0.5),
     data=NULL,  daysAhead=0, Nsim, runs, nthin = runs, proposalOffset) {

     days  <-  lengthOfEpidemic$censor[which(lengthOfEpidemic$censor != "NA")][1]


   
     if(is.null(data) & is.null(attributes(postSample)$data)) {
          warning("need one of either data or augmented data as part of postSample")
     }

	 if(is.null(data)) {
		 data = attributes(postSample)$data[[Dsample]] 
 	}
	 
     Npost = dim(postSample)[1]

     unobservedCases = totalCases = array(NA, c(samples=Npost, days=days + daysAhead, type=3), 
          dimnames=list(NULL,NULL, c("M","S","D")))
    
     for(Dsample in 1:Npost) {
		paramsNew = vecToListParams(postSample[Dsample, ])
		
          dataAug = dataAugment(data, paramsNew) # input arguments are data (data that needs to be augmented) and params
      
          if(!all(c("M","S","D") %in% names(infectParams))) {      # need to double check!
               for(Dtype in c("M","S","D") ) {
                    infectParams[Dtype] =  mean(dataAug[,paste("prob", Dtype, sep="")])
               }
          }
   # day 1, flat prior on N
   
#   prob(N | Y) = prob(Y | N) prob(N) / prob(Y)
# sim cases on the first day
    Dday=1

    medDay = dataAug$med
    infectDay = medDay+dataAug$infect + dataAug$onset
    
    removedDay = dataAug$removed 
    theHosp = dataAug$type != "M"
    removedDay[theHosp]= dataAug[theHosp, "hospital"]
    removedDay = removedDay+ medDay
	
	
	newCases = NULL

 probObs = probSumWeibulls(postSample[Dsample,], x = days - Dday + 0.5, Nsim)
 
 	for(Dtype in c("M", "S", "D")) {

               YobsToday = sum(infectDay == Dday & medDay < days & dataAug$type == Dtype)
               
    		   
			   
               totalCases[Dsample,1,Dtype] = updateBinomialPoisson(priorMean = 10, YobsToday, 
					   probObs$probs[Dtype], 
                    N = YobsToday+10, runs, nthin, proposalOffset)

               unobservedCases[Dsample,1,Dtype] = totalCases[Dsample,1,Dtype]  - YobsToday
   
	newCaseTimes = probObs$sample[[Dtype]][sample(
					dim(probObs$sample[[Dtype]])[2], 
					unobservedCases[Dsample, Dday, Dtype],
					replace=T),
    		]		   
			   
	NAvec = rep(NA, unobservedCases[Dsample, Dday, Dtype] )		
			
   newCasesToday=data.frame(
		   med= rep(Dday, unobservedCases[Dsample, Dday, Dtype]) + newCaseTimes[,"sum"],
 	  	onset =  - newCaseTimes[,"OnsMed"],
    	infect = - newCasesToday[,"InfOns"],
		hospital = NAvec, removed = NAvec, 
		censor =rep( days,unobservedCases[Dsample, Dday, Dtype] ) , 
		type = rep(Dtype,unobservedCases[Dsample, Dday, Dtype] ) , observedType = NAvec, 
		lost=rep(F,unobservedCases[Dsample, Dday, Dtype] )
	)


	newCases = rbind(newCases, 
			newCasesToday[,names(newCases)]
	)
  } # end loop through types
   
# fill in missing data for unobserved cases
  
  newCases = dataAugment(newCases, paramsNew)
  # add new cases to the dataset
  dataAug = rbind(dataAug, newCases[,names(dataAug)])
  
     
          for(Dday in 2:days) { 
    medDay = dataAug$med
    infectDay = medDay+dataAug$infect + dataAug$onset
    
    removedDay = dataAug$removed 
    theHosp = dataAug$type != "M"
    removedDay[theHosp]= dataAug[theHosp, "hospital"]
    removedDay = removedDay+ medDay
    
	NinfectiveToday = sum(infectDay < Dday & removedDay > Dday) 
    
    # get prob of med on days given infected on Dday
    probObs = probSumWeibulls(postSample, x = days - Dday + 0.5, Nsim)  # this will be a vector of length 3, 
      
          # find total cases infected on Dday and had medical visit before days (from dataAug)
  newCases = NULL
      for(Dtype in names(probObs$prob) ) {
                    YobsToday = sum(infectDay == Dday & medDay < days & dataAug$type == Dtype)

                    priorMean = (NinfectiveToday*infectParams["rate"] + infectParams["immigration"]) * infectParams[Dtype]

                    totalCases[Dsample, Dday, Dtype] = 
							updateBinomialPoisson(priorMean = priorMean, YobsToday, probObs$probs[Dtype], 
                         	N = totalCases[Dsample, Dday-1, Dtype], runs, nthin, proposalOffset)
      
                    unobservedCases[Dsample, Dday, Dtype] = totalCases[Dsample, Dday, Dtype] - YobsToday
            
					newCaseTimes = probObs$sample[[Dtype]][sample(
									dim(probObs$sample[[Dtype]])[2], 
									unobservedCases[Dsample, Dday, Dtype],
									replace=T),
    						]		   
			    	
					NAvec = rep(NA, unobservedCases[Dsample, Dday, Dtype] )		
					
    				newCasesToday=data.frame(
		    				med= rep(Dday, unobservedCases[Dsample, Dday, Dtype]) + newCaseTimes[,"sum"],
 	  						onset =  - newCaseTimes[,"OnsMed"],
    						infect = - newCasesToday[,"InfOns"],
					hospital = NAvec, removed = NAvec, 
					censor =rep( days,unobservedCases[Dsample, Dday, Dtype] ) , 
					type = rep(Dtype,unobservedCases[Dsample, Dday, Dtype] ) , observedType = NAvec, 
					lost=rep(F,unobservedCases[Dsample, Dday, Dtype] )
				)
			
			
			newCases = rbind(newCases,	newCasesToday)
			     
     } # end loop through type

# fill in missing data for unobserved cases

newCases = dataAugment(newCases, paramsNew)

dataAugNew = rbind(dataAugNew, newCases[,names(dataAug)])

} # end loop through days
      
          for(Dday in seq(days + 1, length=daysAhead)) {

    		  medDay = dataAug$med
    		  infectDay = medDay+dataAug$infect + dataAug$onset
    		  
    		  removedDay = dataAug$removed 
    		  theHosp = dataAug$type != "M"
    		  removedDay[theHosp]= dataAug[theHosp, "hospital"]
    		  removedDay = removedDay+ medDay
    		  
			  NinfectiveToday = sum(infectDay < Dday & removedDay > Dday) 
              newCases = NULL
               for (Dtype in names(probObs$probs)) {
      
                    totalCases[Dsample, Dday, Dtype] = rpois(1,
                         (NinfectiveToday*infectParams["rate"] + infectParams["immigration"])*infectParams[Dtype])
            
                    unobservedCases[Dsample, Dday, Dtype] = totalCases[Dsample, Dday, Dtype]
            
            
					NAvec = rep(NA, unobservedCases[Dsample, Dday, Dtype] )		
                	
           newCases = rbind(newCases, data.frame(
    infect= rep(Dday, unobservedCases[Dsample, Dday, Dtype]),
        onset = NAvec, med = NAvec, hospital = NAvec, removed = NAvec, 
       censor =NAvec, 
        type = rep(Dtype,unobservedCases[Dsample, Dday, Dtype] ) , observedType = NAvec, 
         lost=rep(F,unobservedCases[Dsample, Dday, Dtype] ) )
		)
           
               } #end loop through type
  		   newCases = dataAugment(newCases, paramsNew)
      	   dataAugNew = rbind(dataAug, newCases) 
			   
      
          } # end loop through days ahead.
      
  
	return(list(unobserved=unobservedCases, total=totalCases)
	)      
}


   
   
plotUnobservedCases <- function(unobservedData, byType) {
  
     unobservedCases = unobservedData
         
     unobservedCasesMild <- t(unobservedCases[,,"M"])
     unobservedCasesSerious <- t(unobservedCases[,,"S"])
     unobservedCasesDeadly <- t(unobservedCases[,,"D"])
         
     Npost = dim(unobservedCasesMild)[2]       
     totalDays = dim(unobservedCasesMild)[1]    
         
     unobservedCasesTotal <- unobservedCasesMild + unobservedCasesSerious + unobservedCasesDeadly
     unobservedCasesTotal <- cbind(unobservedCasesTotal, "dayOfEpidemic" = 1:totalDays)
                 
     unobservedSummaryByType <- rbind(t(unobservedCases[,,"M"]), t(unobservedCases[,,"S"]), t(unobservedCases[,,"D"]))
     # unobservedSummaryByType <- data.frame(unobservedSummaryByType, "epidemic" = rep(c(rep("epidemic", each = days), rep("ahead", each = daysAhead)), 3), "dayOfEpidemic" = rep(1:daysTotal, 3))
     unobservedSummaryByType <- data.frame(unobservedSummaryByType, "infectionType" = c(rep("M", each = totalDays),  rep("S", each = totalDays), rep("D", each = totalDays)))
                            
     for (i in 1:Npost)
          {names(unobservedSummaryByType)[i] <- c(paste("postSample", i, sep = ""))}
        
     if (byType == "TRUE") {                   
               
          matplot(1:totalDays, unobservedCasesMild, type = "l", lty = 1, col = "grey", xlab = "day of epidemic", ylab = "number of cases", pch = 1, ylim = c(0, max(unobservedCasesMild, unobservedCasesSerious, unobservedCasesDeadly)))
                     
          lines(1:totalDays, apply(unobservedCasesMild, 1, mean), col = "red", lwd = 2)
          lines(1:totalDays, apply(unobservedCasesMild, 1, function(x) quantile(x, probs = 0.025)), col = "red", lwd = 2, lty = 2)
          lines(1:totalDays, apply(unobservedCasesMild, 1, function(x) quantile(x, probs = 0.975)), col = "red", lwd = 2, lty = 2)
          
          matplot(1:totalDays, unobservedCasesSerious, type = "l", lty = 1, col = "blue", xlab = "day of epidemic", ylab = "number of cases", pch = 1, add = TRUE)
          
          lines(1:totalDays, apply(unobservedCasesSerious, 1, mean), col = "red", lwd = 2)
          lines(1:totalDays, apply(unobservedCasesSerious, 1, function(x) quantile(x, probs = 0.025)), col = "red", lwd = 2, lty = 2)
          lines(1:totalDays, apply(unobservedCasesSerious, 1, function(x) quantile(x, probs = 0.975)), col = "red", lwd = 2, lty = 2)
         
          matplot(1:totalDays, unobservedCasesDeadly, type = "l", lty = 1, col = "black", xlab = "day of epidemic", ylab = "number of cases", pch = 1, add = TRUE)
          
          lines(1:totalDays, apply(unobservedCasesDeadly, 1, mean), col = "red", lwd = 2)
          lines(1:totalDays, apply(unobservedCasesDeadly, 1, function(x) quantile(x, probs = 0.025)), col = "red", lwd = 1.5, lty = 2)
          lines(1:totalDays, apply(unobservedCasesDeadly, 1, function(x) quantile(x, probs = 0.975)), col = "red", lwd = 1.5, lty = 2)
          
          legend("topright", legend = c("Mild", "Serious", "Deadly"), col = c("grey", "blue", "black"), lty = c(1,1,1))
         
     } else {

          matplot(unobservedCasesTotal[,"dayOfEpidemic"], unobservedCasesTotal[,1:Npost], type = "l", pch = 1, lty = 1, col = "grey", xlab = "day of epidemic", ylab = "number of cases")
     
          lines(unobservedCasesTotal[,"dayOfEpidemic"], apply(unobservedCasesTotal[,1:Npost], 1, mean), col = "black", lwd = 3)
          lines(unobservedCasesTotal[,"dayOfEpidemic"], apply(unobservedCasesTotal[,1:Npost], 1, function(x) quantile(x, probs = 0.025)), col = "black", lwd = 3, lty = 2)
          lines(unobservedCasesTotal[,"dayOfEpidemic"], apply(unobservedCasesTotal[,1:Npost], 1, function(x) quantile(x, probs = 0.975)), col = "black", lwd = 3, lty = 2)
     }

}

