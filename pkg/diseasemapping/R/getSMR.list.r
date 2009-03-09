getSMR.list(popdata, model, casedata = NULL, regionCode = "CSDUID",regionCodeCases = "CSD2006", years = NULL, year.range = NULL,
    area = FALSE, area.scale = 1, ...) {
#  lennon's stuff

    
        isSP = (class(popdata[[1]]) == "SpatialPolygonsDataFrame")
   
    if (is.null(cyears)) {
        cyears = as.integer(names(popdata))
    }        
    if (area & isSP) {
       
            areas <- lapply(lapply(popdata, area), as.numeric)
            for (i in 1:length(popdata)) {
                popdata[[i]]$sqk <- areas[[i]] * area.scale
            }
    }
    
    poplong <- formatPopulation.list(popdata, breaks = attributes(model)$breaks$breaks,
      years = model$xlevels$YEAR, mustAggregate = FALSE,years=years,year.range=year.range)

       
         
         
    ##list if df
    listpop<-tapply(poplong,poplong$YEAR,getSMR,model, casedata = NULL, regionCode =regionCode,
                     regionCodeCases = regionCodeCases, years = NULL, year.range = NULL,
                     area = FALSE, area.scale = 1)

         #if input is  list of sp, return list of sp
        if (isSP) {
            for (i in 1:length(cyears)) {
                popdata[[i]]@data = listpop[[i]]
            }
            listpop<-popdata
         }   
        #else return list of df
  
    listpop
}

