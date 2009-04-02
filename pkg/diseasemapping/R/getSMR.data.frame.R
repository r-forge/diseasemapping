`getSMR` <- function(popdata, model, casedata, regionCode = "CSDUID",
    regionCodeCases = "CSD2006", area = FALSE, area.scale = 1,formatPop=TRUE, ...){
       UseMethod("getSMR")
 }
 
getSMR.data.frame <- function(popdata, model, casedata, regionCode = "CSDUID",
    regionCodeCases = "CSD2006", area = FALSE, area.scale = 1, formatPop=TRUE,...){
#  getSMR(popdata@data, ...)
    if(formatPop){
    poplong <- formatPopulation(popdata, breaks=attributes(model)$breaks$breaks, 
      mustAggregate = FALSE)
     }else{poplong<-popdata}
     # get rid of zero populations,because they dont lead to rates of exactly zero
     poplong=poplong[poplong[,
      grep("^population$", names(poplong), value=T, ignore.case=T)]>0, ]     
    #changes poplong names to be consistent with model
    agevar<-grep("^age$",names(model$xlevels),value=TRUE,ignore.case=TRUE)
    sexvar<-grep("^sex$",names(model$xlevels),value=TRUE,ignore.case=TRUE)
    yearvar<-grep("^year$",names(model$xlevels),value=TRUE,ignore.case=TRUE)

    agevar1<-grep("^age$",names(poplong),value=TRUE,ignore.case=TRUE)
    sexvar1<-grep("^sex$",names(poplong),value=TRUE,ignore.case=TRUE)
    yearvar1<-grep("^year$",names(poplong),value=TRUE,ignore.case=TRUE)
    
    if(length(agevar) & length(agevar1)){
     names(poplong[[agevar1]])=agevar
    }
    if(length(sexvar) & length(sexvar1)){
     names(poplong[[sexvar1]])=sexvar
    }
    if(length(yearvar) & length(yearvar1)){
     names(poplong[[yearvar1]])=yearvar
    }

    
    if (length(model$sexSubset) == 1) {
        poplong = poplong[poplong$sex == model$sexSubset, ]
    }
    
      
      offsetvar<- grep("logpop",names(model$data) ,value=TRUE,ignore.case=TRUE)
    poplong[[offsetvar]] = log(poplong$POPULATION)
#    poplong[[sexvar]]= factor(poplong[[sexvar]])
#    poplong[[agevar]] = factor(poplong[[agevar]])
    #names(poplong) <- tolower(names(poplong))
    for (Dlevel in names(model$xlevels)) {
        alllevels = levels(poplong[[Dlevel]])
        if (!all(alllevels %in% model$xlevels[[Dlevel]])) {
            tokeep = poplong[[Dlevel]] %in% model$xlevels[[Dlevel]]
            poplong = poplong[tokeep, ]
        }
    }
    interactNA <- names(model$coefficients)[is.na(model$coefficients)]
    if (length(interactNA) > 0) {
        interact <- grep(":", interactNA, value = TRUE)
        poplong$param <- paste(paste("age", poplong$age, sep = ""),
            paste("sex", poplong$sex, sep = ""), sep = ":")
        poplong = poplong[!poplong$param %in% interact, ]
        poplong$param <- NULL
    }
    if(length(yearvar)) {   
           agg<-c(yearvar,agevar, sexvar,offsetvar)
           }else{
                agg<-c(agevar, sexvar, offsetvar)
           }

    poplong$expected <- predict(model, poplong[, agg], type = "response")

    
    poplong <- aggregate(poplong$expected, list(poplong[[regionCode]]), sum)
    names(poplong) <- c(regionCode, "expected")

     # remove 'expected' column from population data before merging
     popdata = popdata[,-grep("^observed|^expected", names(popdata))]

     # merge results back into original dataset
     populationData = merge(popdata, poplong, by=regionCode, all.x=T)
     populationData$expected[populationData$expected==0] = NA
    if (area & ("sqk" %in% names(populationData) ) ) {
        populationData$expected_sqk <- 
          populationData$expected/populationData$sqk
    }
     populationData$logExpected = log(populationData$expected)

   
    if (!is.null(casedata)) {
       casedata = formatCases(casedata, ageBreaks=attributes(poplong)$breaks)
       casedata = casedata[
          as.character(casedata[, regionCodeCases]) %in% 
             as.character(poplong[, regionCode]), ]
       casecol = grep("^cases$", names(casedata), value = TRUE,
            ignore.case = TRUE)
       if (!length(casecol)) {
            casecol = "cases"
            cases[, casecol] = 1
       }
       casedata <- aggregate(casedata[[casecol]], 
         list(casedata[[regionCodeCases]]), sum)
       names(casedata) = c(regionCodeCases, "observed")
   
       populationData = merge(populationData, casedata, 
         by.x=regionCode, by.y =regionCodeCases, all.x = TRUE)
        
       populationData$SMR <- populationData$observed/populationData$expected
   }

   populationData
}


