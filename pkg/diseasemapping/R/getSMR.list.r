getSMR.list <- function(popdata, model, casedata = NULL, regionCode = "CSDUID", 
              regionCodeCases = "CSD2006", area = FALSE, area.scale = 1, 
              years = NULL, year.range = NULL, ...){ 
#  lennon's stuff
        #isSP = (class(popdata[[1]]) == "SpatialPolygonsDataFrame")
   
    if (is.null(years)) {
        years = as.integer(names(popdata))
    }        
   # if (area & isSP) {
  #     
  #          areas <- lapply(lapply(popdata, area), as.numeric)
   #         for (i in 1:length(popdata)) {
   #             popdata[[i]]$sqk <- areas[[i]] * area.scale
   #         }
   # }
    
   # poplong <- formatPopulation(popdata, breaks = attributes(model)$breaks$breaks,
    #  years = model$xlevels$YEAR, mustAggregate = FALSE, year.range=year.range)

       
    #ll<-split(poplong,poplong$YEAR)     

  yearVar = NA  # find year var in the model
  caseYearVar = NA
  
  result=list()

  for(Dyear in seq(1, length(popdata))) {
  
    if(ares & class(popdata[[Dyear]]) == "SpatialPolygonsDataFrame") {
      areas <- unlist(lapply(popdata[[Dyear]], area))
      popdata[[Dyear]]$sqk <- areas * area.scale
    }

  
    # scaling factor to convert to person years
    attributes(popdata[[Dyear]])$popScale = 1
    # add year column
    popdata[[Dyear]][,yearVar] = 
       factor(rep(as.character(years[Dyear]), length(popdata[[Dyear]][1])), 
    levels=levels(model$data[,yearVar])) 
    
    caseThisYear = NA # subset just this period's cases, find year col in casedata
    result[[Dyear]] = getSMR(popdata[[Dyear]], caseThisYear,model,
      ...)
    
  }
  names(result) = years
         
  return(result)       
         
    ##list if df
    listpop<-lapply(ll, getSMR, casedata=casedata, model, regionCode =regionCode,
                     regionCodeCases = regionCodeCases, years = years, year.range = year.range,
                     area = area, area.scale = area.scale,formatPop=FALSE)
                                
         #if input is  list of sp, return list of sp
        if (isSP) {
            for (i in 1:length(years)) {
                popdata[[i]]@data = listpop[[i]]
            }
            listpop<-popdata
         }   
        #else return list of df
  
    listpop
}

