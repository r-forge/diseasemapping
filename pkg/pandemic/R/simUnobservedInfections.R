

simUnobservedInfections = function(postSample, lengthOfEpidemic, infectParams = c(rate = 1, immigration = 0.5),
     data=NULL,  daysAhead, Nsim, runs, nthin = runs, proposalOffset) {

     days  <-  lengthOfEpidemic$censor[which(lengthOfEpidemic$censor != "NA")][1]

     paramsNew <- postSample[dim(postSample)[1],]  
      

     paramsNew <- pandemicParams(InfOns = c("mean" = unname(paramsNew["InfOns.mean"]), "shape" = unname(paramsNew["InfOns.shape"]), "scale" = unname(paramsNew["InfOns.zeros"])), 
     OnsMedM = c("mean" = unname(paramsNew["OnsMedM.mean"]), "shape" = unname(paramsNew["OnsMedM.shape"]), "zeros" = unname(paramsNew["OnsMedM.zeros"])), 
     OnsMedS = c("mean" = unname(paramsNew["OnsMedS.mean"]), "shape" = unname(paramsNew["OnsMedS.shape"]), "zeros" = unname(paramsNew["OnsMedS.zeros"])), 
     OnsMedD = c("mean" = unname(paramsNew["OnsMedD.mean"]), "shape" = unname(paramsNew["OnsMedD.shape"]), "zeros" = unname(paramsNew["OnsMedD.zeros"])), 
     MedRec = c("mean" = unname(paramsNew["MedRec.mean"]), "shape" = unname(paramsNew["MedRec.shape"]), "zeros" = unname(paramsNew["MedRec.zeros"]), "lost" = unname(paramsNew["MedRec.lost"])), 
     MedHospS = c("mean" = unname(paramsNew["MedHospS.mean"]), "shape" = unname(paramsNew["MedHospS.shape"]), "zeros" = unname(paramsNew["MedHospS.zeros"])), 
     MedHospD = c("mean" = unname(paramsNew["MedHospD.mean"]), "shape" = unname(paramsNew["MedHospD.shape"]), "zeros" = unname(paramsNew["MedHospD.zeros"])), 
     HospRec = c("mean" = unname(paramsNew["HospRec.mean"]), "shape" = unname(paramsNew["HospRec.shape"]), "zeros" = unname(paramsNew["HospRec.zeros"])), 
     HospDeath = c("mean" = unname(paramsNew["HospDeath.mean"]), "shape" = unname(paramsNew["HospDeath.shape"]), "zeros" = unname(paramsNew["HospDeath.zeros"])), 
     probs = c("M" = unname(paramsNew["probs.M"]), "S" = unname(paramsNew["probs.S"]), "D" = unname(paramsNew["probs.D"])))
      
     if(is.null(data) & is.null(attributes(postSample)$data)) {
          warning("need one of either data or augmented data as part of postSample")
     }

     Npost = dim(postSample)[1]

     unobservedCases = totalCases = array(NA, c(samples=Npost, days=days + daysAhead, type=3), 
          dimnames=list(NULL,NULL, c("M","S","D")))
    
     for(Dsample in 1:Npost) {
          if(is.null(attributes(postSample)$data)) {
               dataAug = dataAugment(data, paramsNew) # input arguments are data (data that needs to be augmented) and params
          } else {
               dataAug = attributes(postSample)$data[[Dsample]]
          }
      
          infectDay = (-1)*dataAug$infect + dataAug$onset
          medDay = dataAug$med
          removedDay = dataAug$removed 
          theHosp = dataAug$type != "M"
          removedDay[theHosp]= dataAug[theHosp, "hospital"]
          removedDay = removedDay+ medDay

          if(!all(c("M","S","D") %in% names(infectParams))) {      # need to double check!
               for(Dtype in c("M","S","D") ) {
                    infectParams[Dtype] =  mean(dataAug[,paste("prob", Dtype, sep="")])
               }
          }
   # day 1, flat prior on N
   
#   prob(N | Y) = prob(Y | N) prob(N) / prob(Y)
          Dday=1
 newCases = NULL
          for(Dtype in c("M", "S", "D")) {

               YobsToday = sum(infectDay == Dday & medDay < days & dataAug$type == Dtype)
               
               probObs = probSumWeibulls(postSample, x = days - Dday + 0.5, Nsim)$prob  
               totalCases[Dsample,1,Dtype] = updateBinomialPoisson(priorMean = 10, YobsToday, probObs[Dtype], 
                    N = YobsToday+10, runs, nthin, proposalOffset)

               unobservedCases[Dsample,1,Dtype] = totalCases[Dsample,1,Dtype]  - YobsToday
   
   newCases = rbind(newCases, data.frame(
    infect= rep(Dday, unobservedCases[Dsample, Dday, Dtype]),
         onset = NA, med = NA, hospital = NA, removed = NA, 
         censor =rep( days,unobservedCases[Dsample, Dday, Dtype] ) , 
         type = rep(Dtype,unobservedCases[Dsample, Dday, Dtype] ) , observedType = NA, 
          lost=rep(F,unobservedCases[Dsample, Dday, Dtype] ) )
      )
          }
   
 
   newCases = dataAugment(newCases, paramsNew)
  dataAugNew = rbind(dataAug, newCases)

     
          for(Dday in 2:days) { 
          # get prob of med on days given infected on Dday
               probObs = probSumWeibulls(postSample, x = days - Dday + 0.5, Nsim)  # this will be a vector of length 3, 
      
          # find total cases infected on Dday and had medical visit before days (from dataAug)
  
                                   
               NinfectiveToday = sum(infectDay < Dday & removedDay > Dday) 
 
               for(Dtype in names(probObs) ) {
                    YobsToday = sum(infectDay == Dday & medDay < days & dataAug$type == Dtype)

                    priorMean = (NinfectiveToday*infectParams["rate"] + infectParams["immigration"]) * infectParams[Dtype]

                    totalCases[Dsample, Dday, Dtype] = updateBinomialPoisson(priorMean = priorMean, YobsToday, probObs[Dtype], 
                         N = totalCases[Dsample, Dday-1, Dtype], runs, nthin, proposalOffset)
      
                    unobservedCases[Dsample, Dday, Dtype] = totalCases[Dsample, Dday, Dtype] - YobsToday
            
              newCases = rbind(newCases, data.frame(
      infect= rep(Dday, unobservedCases[Dsample, Dday, Dtype]),
           onset = NA, med = NA, hospital = NA, removed = NA, 
           censor =rep( days,unobservedCases[Dsample, Dday, Dtype] ) , 
           type = rep(Dtype,unobservedCases[Dsample, Dday, Dtype] ) , observedType = NA, 
            lost=rep(F,unobservedCases[Dsample, Dday, Dtype] ) )
          
               }
     # loop through days ahead
     
   newCases = data.frame(infect= rep(Dday, unobservedCases[Dsample, Dday, Dtype]),
          onset = NA, med = NA, hospital = NA, removed = NA, censor = NA, type = Dtype, observedType = NA)
       newCases = dataAugment(newCases, paramsNew)
      dataAugNew = rbind(dataAug, newCases)
     

          }
        
          newCases = dataAugment(newCases, paramsNew)
       dataAugNew = rbind(dataAug, newCases) 
      
          for(Dday in seq(days + 1, length=daysAhead)) {

               NinfectiveToday = sum(infectDay < Dday & removedDay > Dday) 

               for (Dtype in names(probObs)) {
      
                    totalCases[Dsample, Dday, Dtype] = rpois(1,
                         (NinfectiveToday*infectParams["rate"] + infectParams["immigration"])*infectParams[Dtype])
            
                    unobservedCases[Dsample, Dday, Dtype] = totalCases[Dsample, Dday, Dtype]
            
            
              
           newCases = rbind(newCases, data.frame(
    infect= rep(Dday, unobservedCases[Dsample, Dday, Dtype]),
        onset = NA, med = NA, hospital = NA, removed = NA, 
       censor =rep( days,unobservedCases[Dsample, Dday, Dtype] ) , 
        type = rep(Dtype,unobservedCases[Dsample, Dday, Dtype] ) , observedType = NA, 
         lost=rep(F,unobservedCases[Dsample, Dday, Dtype] ) )

           
               }

      
          }
      
      
   newCases = data.frame(infect= rep(Dday, unobservedCases[Dsample, Dday, Dtype]),
          onset = NA, med = NA, hospital = NA, removed = NA, censor = NA, type = Dtype, observedType = NA)
       newCases = dataAugment(newCases, paramsNew)
      dataAugNew = rbind(dataAugNew, newCases)
       
       }
   
  newCases = dataAugment(newCases, paramsNew)
      dataAugNew = rbind(dataAug, newCases) 


  
          unobservedCases
      
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

