getBreaks <- function(colNames) {
popColumns <- grep("^(m|f|male|female)[[:digit:]]+(_|-|plus|\\+)[[:digit:]]*$", colNames, value=T, ignore.case=T)
ageGroups <- gsub("^(m|f|male|female)", "", popColumns, ignore.case=TRUE)
ageGroups <- gsub("(\\+|plus)", "_Inf", ageGroups, ignore.case=TRUE)
ageLower <- as.numeric(gsub("(_|-)([[:digit:]]+|Inf)$", "", ageGroups) )
ageUpper <- as.numeric(gsub("^[[:digit:]]+(_|-)", "", ageGroups) )
currentbreaks <- c(sort(unique(ageLower)), max(ageUpper))

sex <- substr(popColumns,1,1)
ageCut <- as.character(cut(ageLower, breaks=currentbreaks, right=F)) # I
return(list(breaks=currentbreaks, names=ageCut, age=ageLower, sex=sex, 
	oldNames=popColumns, 
	newNames=paste(sex, ageLower, sep=".")))
}