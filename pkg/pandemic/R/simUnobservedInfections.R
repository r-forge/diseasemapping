


# description: this function runs a Binomial-Poisson MCMC to calculate the number of unobserved infections and 
     # forecases the number of unobserved infections beyond the simulated epidemic

# input: 
     # postSample: output from mcmcPandemic
     # lengthOfEpidemic: output from simEpidemic
     # infectParams: a vector with the rate and immigration parameters 
     # data: NULL (if getting data from mcmcPandemic) or enter a data set to run through dataAugment function
     # daysAhead: number of days beyond the end of the simulated epidemic that you want to forecast
     # Nsim: number of random variables simulated to calculate the probability in the probSumWeibulls function
     # runs: number of runs of the Binomial Poisson MCMC
     # nthin: thinning parameter (when runs = nthin, we only retain the last run)
     # proposalOffset: proposal offset for the MCMC simulation

simUnobservedInfections = function(postSample, lengthOfEpidemic, infectParams = c(rate = 1, immigration = 0.5), data = NULL,
  daysAhead = 2, Nsim = 1000, runs = 100, nthin = 100, proposalOffset = 0.2) {

     result=list()
     days = lengthOfEpidemic$censor[which(lengthOfEpidemic$censor != "NA")][1]

     if(is.null(data) & is.null(attributes(postSample)$data)) {
          warning("need one of either data or augmented data as part of postSample")
     }

     Npost = dim(postSample)[1]
     unobservedCases = totalCases = array(NA, c(samples=Npost, days=days + daysAhead, type=3),
          dimnames=list(NULL,NULL, c("M","S","D")))

     for(Dsample in 1:Npost) {
     #print(Dsample)

     dataAug = NULL

          if(is.null(data)) {
               dataAug = attributes(postSample)$data[[Dsample]]
          } else {
               dataAug = dataAugment(data, paramsNew) # input arguments are data (data that needs to be augmented) and params
          }

          paramsNew = vecParamsToList(postSample[Dsample, ])

          if(!all(c("M","S","D") %in% names(infectParams))) {

               for(Dtype in c("M","S","D") ) {
                    infectParams[Dtype] =  mean(dataAug[,paste("prob", Dtype, sep="")])
               }

          }
# day 1, flat prior on N
# prob(N | Y) = prob(Y | N) prob(N) / prob(Y)
# sim cases on the first day
          Dday=1
          medDay = dataAug$med
          infectDay = medDay + dataAug$infect
          removedDay = dataAug$removed
          theHosp = dataAug$type != "M"
          removedDay[theHosp]= dataAug[theHosp, "hospital"]
          removedDay = removedDay+ medDay

          #testLength = as.numeric(length(medDay) == length(infectDay) & length(medDay) == length(removedDay))
          #print(testLength)

          newCases = NULL
          probObs = probSumWeibulls(postSample[Dsample,], x = days - Dday + 0.5, Nsim)

          for(Dtype in c("M", "S", "D")) {
               YobsToday = sum(infectDay == Dday & medDay < days & dataAug$type == Dtype)
               totalCases[Dsample,1,Dtype] = updateBinomialPoisson(priorMean = 1, Y = YobsToday,
                    prob = 1 - probObs$prob[Dtype], runs = runs, nthin = nthin,
                    proposalOffset = proposalOffset)
                    unobservedCases[Dsample,1,Dtype] = totalCases[Dsample,1,Dtype]  - YobsToday
                    NAvec = rep(NA, unobservedCases[Dsample, 1, Dtype])

               if (unobservedCases[Dsample, 1, Dtype] > 0) {
                    newCaseTimes = probObs$sample[[Dtype]][sample(dim(probObs$sample[[Dtype]])[1],
                    unobservedCases[Dsample, 1, Dtype], replace=T), ]

                    if (unobservedCases[Dsample, Dday, Dtype] == 1) {
                         newCasesToday=data.frame(med= rep(Dday, unobservedCases[Dsample, 1, Dtype]) +
                              newCaseTimes["sum"], onset = -newCaseTimes["OnsMed"], infect = -newCaseTimes["InfOns"],
                              hospital = NAvec, removed = NAvec, censor =rep(days, unobservedCases[Dsample, 1, Dtype]),
                              type = rep(Dtype,unobservedCases[Dsample, 1, Dtype]), observedType = rep(Dtype,unobservedCases[Dsample, Dday, Dtype]),
                              lost=rep(FALSE,unobservedCases[Dsample, 1, Dtype]), age=NAvec)
                    }

                    if (unobservedCases[Dsample, Dday, Dtype] > 1) {
                         newCasesToday=data.frame(med = rep(Dday, unobservedCases[Dsample, 1, Dtype]) +
                              newCaseTimes[,"sum"], onset = -newCaseTimes[,"OnsMed"], infect = -newCaseTimes[,"InfOns"],
                              hospital = NAvec, removed = NAvec, censor =rep(days,unobservedCases[Dsample, 1, Dtype]),
                              type = rep(Dtype,unobservedCases[Dsample, 1, Dtype]), observedType = rep(Dtype,unobservedCases[Dsample, Dday, Dtype]),
                              lost = rep(FALSE,unobservedCases[Dsample, 1, Dtype]), age=NAvec)
                    }

                    newCases = rbind(newCases, newCasesToday[,names(newCasesToday)])

               } else {

                    newCasesToday = NULL
               }
          } # end loop through types
  # fill in missing data for unobserved cases
          if (!is.null(newCases)) {
               newCases = dataAugment(newCases, paramsNew)
               dataAug = rbind(dataAug, newCases[,names(dataAug)])
          }

          medDay = dataAug$med
          infectDay = medDay + dataAug$infect
          infectType = dataAug$type
          removedDay = dataAug$removed
          theHosp = dataAug$type != "M"
          removedDay[theHosp] = dataAug[theHosp, "hospital"]
          removedDay = removedDay+ medDay

          #testLength1 = as.numeric(length(medDay) == length(infectDay) & length(medDay) == length(removedDay))
          #print(testLength1)

          for(Dday in 2:days) {
          #print(Dday)

               NinfectiveToday = sum(infectDay < Dday & removedDay > Dday)
# get prob of med on days given infected on Dday
               probObs = probSumWeibulls(postSample[Dsample,], x = days - Dday + 0.5, Nsim)
# this will be a vector of length 3,
# find total cases infected on Dday and had medical visit before days (from dataAug)
               newCases = NULL

               for(Dtype in names(probObs$prob)) {
               #print(Dtype)

                    YobsToday = sum(infectDay == Dday & medDay < days & infectType == Dtype)

                    #testLength3 = as.numeric(length(infectDay) == length(medDay) & length(infectDay) == length(infectType))
                    #print(testLength3)

                    priorMean = (NinfectiveToday*infectParams["rate"] + infectParams["immigration"]) * infectParams[Dtype]
                    totalCases[Dsample, Dday, Dtype] = updateBinomialPoisson(priorMean = priorMean, Y=YobsToday,
                         prob=1-probObs$prob[Dtype], runs=runs, nthin=nthin, proposalOffset=proposalOffset)
                    unobservedCases[Dsample, Dday, Dtype] = totalCases[Dsample, Dday, Dtype] - YobsToday


                    if (unobservedCases[Dsample, Dday, Dtype] > 0) {

                         newCaseTimes = (probObs$sample[[Dtype]][sample(dim(probObs$sample[[Dtype]])[1],
                              unobservedCases[Dsample, Dday, Dtype], replace=T), ])
                         NAvec = rep(NA, unobservedCases[Dsample, Dday, Dtype])

                         if (unobservedCases[Dsample, Dday, Dtype] == 1) {

                              newCasesToday=data.frame(med = rep(Dday, unobservedCases[Dsample, Dday, Dtype]) + newCaseTimes["sum"],
                                   onset = -newCaseTimes["OnsMed"], infect = -newCaseTimes["InfOns"], hospital = NAvec, removed = NAvec,
                                   censor = rep(days,unobservedCases[Dsample, Dday, Dtype]), type = rep(Dtype,unobservedCases[Dsample, Dday, Dtype]),
                                   observedType = rep(Dtype,unobservedCases[Dsample, Dday, Dtype]), 
                                   lost=rep(F,unobservedCases[Dsample, Dday, Dtype]), age=NAvec)


                              #newCases = rbind(newCases, newCasesToday)

                         }

                         if (unobservedCases[Dsample, Dday, Dtype] > 1) {

                              newCasesToday = data.frame(med = rep(Dday, unobservedCases[Dsample, Dday, Dtype]) + newCaseTimes[,"sum"],
                                   onset = -newCaseTimes[,"OnsMed"], infect = -newCaseTimes[,"InfOns"], hospital = NAvec, removed = NAvec,
                                   censor = rep(days,unobservedCases[Dsample, Dday, Dtype]), type = rep(Dtype,unobservedCases[Dsample, Dday, Dtype]),
                                   observedType = rep(Dtype,unobservedCases[Dsample, Dday, Dtype]), 
                                   lost=rep(F,unobservedCases[Dsample, Dday, Dtype]), age=NAvec)  

                              #newCases = rbind(newCases, newCasesToday)

                         }

                    } else {

                         newCasesToday = NULL

                    }

                    newCases = rbind(newCases, newCasesToday)

               } # end loop through type

               if (!is.null(newCases)) {

                    newCases = try(dataAugment(newCases, paramsNew))

               }

               dataAug = rbind(dataAug, newCases[,names(dataAug)])

          }  # end loop through days
          #print(dataAug)

          medDay = dataAug$med
          infectDay = medDay + dataAug$infect
          removedDay = dataAug$removed
          theHosp = dataAug$type != "M"
          removedDay[theHosp]= dataAug[theHosp, "hospital"]
          removedDay = removedDay+ medDay

          #testLength2 = as.numeric(length(medDay) == length(infectDay) & length(medDay) == length(removedDay))
          #print(testLength2)

          for(Dday in seq(days + 1, length=daysAhead)) {
          #print(Dday)

               NinfectiveToday = sum(infectDay < Dday & removedDay > Dday)
               newCases = NULL

               for (Dtype in names(probObs$prob)) {
               #print(Dtype)

                    totalCases[Dsample, Dday, Dtype] = rpois(1, (NinfectiveToday*infectParams["rate"] + infectParams["immigration"])*infectParams[Dtype])
                    unobservedCases[Dsample, Dday, Dtype] = totalCases[Dsample, Dday, Dtype]
                    
                      if (unobservedCases[Dsample, Dday, Dtype] > 0) {
                    
               
                    NAvec = rep(NA, unobservedCases[Dsample, Dday, Dtype])
                    newCasesToday = data.frame(infect= rep(Dday, unobservedCases[Dsample, Dday, Dtype]),
                         onset = NAvec, med = NAvec, hospital = NAvec, removed = NAvec, censor =NAvec,
                         type = rep(Dtype,unobservedCases[Dsample, Dday, Dtype]) , observedType = rep(Dtype,unobservedCases[Dsample, Dday, Dtype]),
                         lost = rep(FALSE,unobservedCases[Dsample, Dday, Dtype]), age = NAvec)
                    #newCases = rbind(newCases,	newCasesToday)
         
                    
             

                    } else {
                    
                      newCasesToday = NULL
                    
                    
                    }

                    newCases = rbind(newCases, newCasesToday)

               } #end loop through type

                if (!is.null(newCases)) {
                    newCases = try(dataAugment(newCases, paramsNew))
                    dataAug = rbind(dataAug, newCases[,names(dataAug)])
               } # end loop through days ahead.


               #if (!is.null(newCases)) {
                #    newCases = dataAugment(newCases, paramsNew)
                #    dataAug = rbind(dataAug, newCases)
               #} # end loop through days ahead.

          }

          result[[Dsample]] = dataAug

     } # end loop through number of posterior samples

     list("unobserved" = unobservedCases, "total" = totalCases, "sample" = result)
 }
 

# input:
      # cases: output from simUnobservedInfections
      # byType: TRUE (if you want to plot the unobserved cases by type) and FALSE (if you want a plot of all unobserved cases) 
 
plotUnobservedCases <- function(cases, byType) {
  
     unobservedCases = cases$unobserved
         
     unobservedCasesMild <- t(unobservedCases[,,"M"])
     unobservedCasesSerious <- t(unobservedCases[,,"S"])
     unobservedCasesDeadly <- t(unobservedCases[,,"D"])
         
     Npost = dim(unobservedCasesMild)[2]       
     totalDays = dim(unobservedCasesMild)[1]    
         
     unobservedCasesTotal <- unobservedCasesMild + unobservedCasesSerious + unobservedCasesDeadly
     unobservedCasesTotal <- cbind(unobservedCasesTotal, "dayOfEpidemic" = 1:totalDays)
                 
     unobservedSummaryByType <- rbind(t(unobservedCases[,,"M"]), t(unobservedCases[,,"S"]), t(unobservedCases[,,"D"]))
     unobservedSummaryByType <- data.frame(unobservedSummaryByType, "infectionType" = c(rep("M", each = totalDays),  rep("S", each = totalDays), rep("D", each = totalDays)))
                            
     for (i in 1:Npost)
          {names(unobservedSummaryByType)[i] <- c(paste("postSample", i, sep = ""))}
        
     if (byType == "TRUE") {                   
               
          matplot(1:totalDays, unobservedCasesMild, type = "l", lty = 1, col = "grey", xlab = "day of epidemic", ylab = "number of cases", pch = 1, ylim = c(0, max(unobservedCasesMild, unobservedCasesSerious, unobservedCasesDeadly)))
                     
          lines(1:totalDays, apply(unobservedCasesMild, 1, mean), col = "purple", lwd = 2)
          lines(1:totalDays, apply(unobservedCasesMild, 1, function(x) quantile(x, probs = 0.025)), col = "purple", lwd = 2, lty = 2)
          lines(1:totalDays, apply(unobservedCasesMild, 1, function(x) quantile(x, probs = 0.975)), col = "purple", lwd = 2, lty = 2)
          
          matplot(1:totalDays, unobservedCasesSerious, type = "l", lty = 1, col = "blue", xlab = "day of epidemic", ylab = "number of cases", pch = 1, add = TRUE)
          
          lines(1:totalDays, apply(unobservedCasesSerious, 1, mean), col = "yellow", lwd = 2)
          lines(1:totalDays, apply(unobservedCasesSerious, 1, function(x) quantile(x, probs = 0.025)), col = "yellow", lwd = 2, lty = 2)
          lines(1:totalDays, apply(unobservedCasesSerious, 1, function(x) quantile(x, probs = 0.975)), col = "yellow", lwd = 2, lty = 2)
         
          matplot(1:totalDays, unobservedCasesDeadly, type = "l", lty = 1, col = "black", xlab = "day of epidemic", ylab = "number of cases", pch = 1, add = TRUE)
          
          lines(1:totalDays, apply(unobservedCasesDeadly, 1, mean), col = "pink", lwd = 2)
          lines(1:totalDays, apply(unobservedCasesDeadly, 1, function(x) quantile(x, probs = 0.025)), col = "pink", lwd = 1.5, lty = 2)
          lines(1:totalDays, apply(unobservedCasesDeadly, 1, function(x) quantile(x, probs = 0.975)), col = "pink", lwd = 1.5, lty = 2)
          
          legend("topleft", legend = c("Mild", "Serious", "Deadly"), col = c("grey", "blue", "black"), lty = c(1,1,1))
         
     } else {

          matplot(unobservedCasesTotal[,"dayOfEpidemic"], unobservedCasesTotal[,1:Npost], type = "l", pch = 1, lty = 1, col = "grey", xlab = "day of epidemic", ylab = "number of cases")
     
          lines(unobservedCasesTotal[,"dayOfEpidemic"], apply(unobservedCasesTotal[,1:Npost], 1, mean), col = "black", lwd = 3)
          lines(unobservedCasesTotal[,"dayOfEpidemic"], apply(unobservedCasesTotal[,1:Npost], 1, function(x) quantile(x, probs = 0.025)), col = "black", lwd = 3, lty = 2)
          lines(unobservedCasesTotal[,"dayOfEpidemic"], apply(unobservedCasesTotal[,1:Npost], 1, function(x) quantile(x, probs = 0.975)), col = "black", lwd = 3, lty = 2)
     }

}

# input:
      # cases: output from simUnobservedInfections
      # byType: TRUE (if you want to plot the total cases by type) and FALSE (if you want a plot of all cases) 

plotTotalCases <- function(cases, byType) {
  
     totalCases = cases$total
         
     totalCasesMild <- t(totalCases[,,"M"])
     totalCasesSerious <- t(totalCases[,,"S"])
     totalCasesDeadly <- t(totalCases[,,"D"])
         
     Npost = dim(totalCasesMild)[2]       
     totalDays = dim(totalCasesMild)[1]    
         
     totalCasesTotal <- totalCasesMild + totalCasesSerious + totalCasesDeadly
     totalCasesTotal <- cbind(totalCasesTotal, "dayOfEpidemic" = 1:totalDays)
                 
     totalSummaryByType <- rbind(t(totalCases[,,"M"]), t(totalCases[,,"S"]), t(totalCases[,,"D"]))
     totalSummaryByType <- data.frame(totalSummaryByType, "infectionType" = c(rep("M", each = totalDays),  rep("S", each = totalDays), rep("D", each = totalDays)))
                            
     for (i in 1:Npost)
          {names(totalSummaryByType)[i] <- c(paste("postSample", i, sep = ""))}
        
     if (byType == "TRUE") {                   
               
          matplot(1:totalDays, totalCasesMild, type = "l", lty = 1, col = "grey", xlab = "day of epidemic", ylab = "number of cases", pch = 1, ylim = c(0, max(totalCasesMild, totalCasesSerious, totalCasesDeadly)))
                     
          lines(1:totalDays, apply(totalCasesMild, 1, mean), col = "purple", lwd = 2)
          lines(1:totalDays, apply(totalCasesMild, 1, function(x) quantile(x, probs = 0.025)), col = "purple", lwd = 2, lty = 2)
          lines(1:totalDays, apply(totalCasesMild, 1, function(x) quantile(x, probs = 0.975)), col = "purple", lwd = 2, lty = 2)
          
          matplot(1:totalDays, totalCasesSerious, type = "l", lty = 1, col = "blue", xlab = "day of epidemic", ylab = "number of cases", pch = 1, add = TRUE)
          
          lines(1:totalDays, apply(totalCasesSerious, 1, mean), col = "yellow", lwd = 2)
          lines(1:totalDays, apply(totalCasesSerious, 1, function(x) quantile(x, probs = 0.025)), col = "yellow", lwd = 2, lty = 2)
          lines(1:totalDays, apply(totalCasesSerious, 1, function(x) quantile(x, probs = 0.975)), col = "yellow", lwd = 2, lty = 2)
         
          matplot(1:totalDays, totalCasesDeadly, type = "l", lty = 1, col = "black", xlab = "day of epidemic", ylab = "number of cases", pch = 1, add = TRUE)
          
          lines(1:totalDays, apply(totalCasesDeadly, 1, mean), col = "pink", lwd = 2)
          lines(1:totalDays, apply(totalCasesDeadly, 1, function(x) quantile(x, probs = 0.025)), col = "pink", lwd = 1.5, lty = 2)
          lines(1:totalDays, apply(totalCasesDeadly, 1, function(x) quantile(x, probs = 0.975)), col = "pink", lwd = 1.5, lty = 2)
          
          legend("topleft", legend = c("Mild", "Serious", "Deadly"), col = c("grey", "blue", "black"), lty = c(1,1,1))
         
     } else {

          matplot(totalCasesTotal[,"dayOfEpidemic"], totalCasesTotal[,1:Npost], type = "l", pch = 1, lty = 1, col = "grey", xlab = "day of epidemic", ylab = "number of cases")
     
          lines(totalCasesTotal[,"dayOfEpidemic"], apply(totalCasesTotal[,1:Npost], 1, mean), col = "black", lwd = 3)
          lines(totalCasesTotal[,"dayOfEpidemic"], apply(totalCasesTotal[,1:Npost], 1, function(x) quantile(x, probs = 0.025)), col = "black", lwd = 3, lty = 2)
          lines(totalCasesTotal[,"dayOfEpidemic"], apply(totalCasesTotal[,1:Npost], 1, function(x) quantile(x, probs = 0.975)), col = "black", lwd = 3, lty = 2)
     }

}