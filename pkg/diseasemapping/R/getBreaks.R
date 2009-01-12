function (colNames, breaks=NULL) 
{
    if (!length(breaks)){
    popColumns <- grep("^(m|f|male|female)[[:digit:]]+(_|-|plus|\\+)[[:digit:]]*$", 
        colNames, value = T, ignore.case = T)
    ageGroups <- gsub("^(m|f|male|female)", "", popColumns, ignore.case = TRUE)
    ageGroups <- gsub("(\\+|plus)", "_Inf", ageGroups, ignore.case = TRUE)
    ageLower <- as.numeric(gsub("(_|-)([[:digit:]]+|Inf)$", "", 
        ageGroups))
    ageUpper <- as.numeric(gsub("^[[:digit:]]+(_|-)", "", ageGroups))
    currentbreaks <- c(sort(unique(ageLower)), max(ageUpper))
    sex <- substr(popColumns, 1, 1)
    ageCut <- as.character(cut(ageLower, breaks = currentbreaks, 
        right = F))
    ageBreaks = list(breaks = currentbreaks, names = ageCut, age = ageLower, 
        sex = sex, oldNames = popColumns, newNames = paste(sex, 
            ageLower, sep = "."))
}else(
    poplong = reshape(popdata, varying = ageBreaks$oldNames, 
        direction = "long", v.names = "POPULATION", timevar = "GROUP", 
        times = ageBreaks$newNames)
    a <- numeric()
    age <- substr(poplong$GROUP, 3, 4)
    for (i in 1:length(age)) {
        a[i] <- which(age[i] == unique(ageBreaks$age))
        poplong$cutAge[i] <- ageBreaks$names[a[i]]
    }
    popa <- aggregate(poplong$POPULATION, list(cutAge=poplong$cutAge), sum)
    names(popa)[names(popa) == "x"] = "POPULATION"
}
