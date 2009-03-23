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
    #oldnames<-c(agevar1,sexvar1,yearvar1)
#    newnames<-c(agevar,sexvar,yearvar)
#
#    #names(newnames) = oldnames
#
#    tochange=which(names(poplong) %in% oldnames)
#    names(poplong)[tochange]<-newnames[names(poplong)[tochange] ]
    
    if (length(model$sexSubset) == 1) {
        poplong = poplong[poplong$sex == model$sexSubset, ]
    }
    
      poplong <- poplong[poplong$POPULATION > 0, ]
      
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

    if (area) {
        poplong$expected_sqk <- poplong$expected/poplong$sqk
    }
    
   
    if (!is.null(casedata)) {
       casedata = formatCases(casedata, ageBreaks=attributes(poplong)$breaks)
        casedata = casedata[as.character(casedata[, regionCodeCases]) %in% as.character(poplong[, regionCode]), ]
        casecol = grep("^cases$", names(casedata), value = TRUE,
            ignore.case = TRUE)
        if (!length(casecol)) {
            casecol = "cases"
            cases[, casecol] = 1
        }
    casedata <- aggregate(casedata[[casecol]], list(casedata[[regionCodeCases]]), sum)
    names(casedata) = c(regionCodeCases, "observed")
        populationData = merge(poplong, casedata, by.x=regionCode, by.y =regionCodeCases, all.x = TRUE)
        populationData$expected[is.na(populationData$observed)] = 0  
        populationData$logExpected = log(populationData$expected)
        if (!is.null(populationData$observed)) {
            populationData$SMR <- populationData$observed/populationData$expected
        }
    }

populationData
    }


