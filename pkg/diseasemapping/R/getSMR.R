#return a list of dataframes if multiple census are given
#Need case file has region in it

`getSMR` <-
function (model, popdata, casedata = NULL, regionCode = "CSDUID", 
    regionCodeCases = "CSD2006", cyears = NULL, year.range = NULL, 
    area = FALSE, area.scale = 1) 
{
    morethanoneyear = class(popdata) == "list"
    if (morethanoneyear) {
        isSP = (class(popdata[[1]]) == "SpatialPolygonsDataFrame")
    }
    else {
        isSP = (class(popdata) == "SpatialPolygonsDataFrame")
    }
    if (is.null(cyears) & morethanoneyear) {
        cyears = as.integer(names(popdata))
    }
    if (area & isSP) {
        if (morethanoneyear) {
            areas <- lapply(lapply(popdata, area), as.numeric)
            for (i in 1:length(popdata)) {
                popdata[[i]]$sqk <- areas[[i]] * area.scale
            }
        }
        else {
            popdata$sqk <- c(unlist(area(popdata))) * area.scale
        }
    }
    poplong <- formatPopulation(popdata, breaks = model$breaks, years = model$years)
    if (is.null(year.range) & morethanoneyear) {
        year.range = range(poplong$YEAR)
    }
    if (length(model$sexSubset) == 1) {
        poplong = poplong[poplong$sex == model$sexSubset, ]
    }
    if (morethanoneyear) {
        times <- c(year.range[1], cyears, year.range[2])
        inter <- diff(times)/2
        nseq <- 1:length(inter) - 1
        mseq <- 2:length(inter)
        interval <- inter[mseq] + inter[nseq]
        names(interval) <- names(table(poplong$YEAR))
        poplong$yearsForCensus = interval[as.character(poplong$YEAR)]
        poplong$POPULATION = poplong$POPULATION * poplong$yearsForCensus
        poplong$YEAR = factor(poplong$YEAR, levels = unique(poplong$YEAR))
        poplong$YEAR = factor(poplong$YEAR)
    }
    poplong <- poplong[poplong$POPULATION > 0, ]
    poplong$LOGPOP = log(poplong$POPULATION)
    poplong$SEX = factor(poplong$SEX)
    poplong$AGE = factor(poplong$AGE)
    names(poplong) <- tolower(names(poplong))
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
    if (morethanoneyear) {
        agg <- c("year", "age", "sex", "logpop")
    }
    else {
        agg <- c("age", "sex", "logpop")
    }
    poplong$expected <- predict(model, poplong[, agg], type = "response")
    if (area) {
        poplong$expected_sqk <- poplong$expected/poplong$sqk
    }
    regionCode <- tolower(regionCode)
    if (morethanoneyear) {
        if (area) {
            newpop <- aggregate(poplong[, c("expected", "expected_sqk")], 
                list(Region = poplong[[regionCode]], Year = poplong$year), 
                sum)
        }
        else {
            newpop <- aggregate(poplong[, "expected"], list(Region = poplong[[regionCode]], 
                Year = poplong$year), sum)
            names(newpop)[names(newpop) == "x"] = "expected"
        }
    }
    else {
        if (area) {
            newpop <- aggregate(poplong[, c("expected", "expected_sqk")], 
                list(Region = poplong[[regionCode]]), sum)
        }
        else {
            newpop <- aggregate(poplong[, "expected"], list(Region = poplong[[regionCode]]), 
                sum)
            names(newpop)[names(newpop) == "x"] = "expected"
        }
    }
    newpop$logExpected = log(newpop$expected)
    populationData <- newpop
    if (!is.null(casedata)) {
        casedata <- formatCases(casedata)
        casedata = casedata[as.character(casedata[, regionCodeCases]) %in% 
            as.character(populationData[, "Region"]), ]
        casecol = grep("^cases$", names(casedata), value = TRUE, 
            ignore.case = TRUE)
        if (!length(casecol)) {
            casecol = "cases"
            cases[, casecol] = 1
        }
        if (morethanoneyear) {
            newcase <- aggregate(casedata[[casecol]], list(Region = casedata[[toupper(regionCodeCases)]], 
                Year = casedata$YEAR), sum)
        }
        else {
            newcase <- aggregate(casedata[[casecol]], list(Region = casedata[[toupper(regionCodeCases)]]), 
                sum)
        }
        names(newcase)[names(newcase) == "x"] = "Observed"
        populationData = merge(newpop, newcase, all.x = TRUE)
        populationData$Observed[is.na(populationData$Observed)] = 0
        if (!is.null(populationData$Observed)) {
            populationData$SMR <- populationData$Observed/populationData$expected
        }
    }
    if (morethanoneyear) {
        listpop <- list()
        for (i in 1:length(cyears)) {
            listpop[[i]] <- merge(popdata[[i]], populationData[populationData$Year == 
                cyears[i], ], by.x = grep(regionCode, names(popdata[[i]]), 
                ignore.case = TRUE, value = TRUE), by.y = "Region", 
                all.x = TRUE)
        }
        if (isSP) {
            for (i in 1:length(cyears)) {
                popdata[[i]]@data = listpop[[i]]
            }
            listpop = popdata
        }
    }
    else {
        listpop = merge(popdata@data, populationData, by.x = grep(paste("^", 
            regionCode, "$", sep = ""), names(popdata@data), 
            ignore.case = TRUE, value = TRUE), by.y = "Region", 
            all.x = TRUE)
        if (isSP) {
            popdata@data = listpop
            listpop = popdata
        }
    }
    listpop
}

