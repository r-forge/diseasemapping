formatPopulation <- function(popdata, aggregate.by=NULL, breaks=NULL, ...) {
  UseMethod("formatPopulation")
}

`formatPopulation.data.frame` <-
function(popdata, aggregate.by=NULL, breaks=NULL, ...) {

#popdata <- popdata@data
ageBreaks = getBreaks(names(popdata), breaks)

####reshape the popdata:
poplong = reshape(popdata,  varying=ageBreaks$oldNames, direction="long",
	v.names="POPULATION", timevar="GROUP", times = ageBreaks$newNames)
# create age and sex variables
agecol = grep("^age$", names(poplong), value=TRUE, ignore.case=TRUE)
sexcol = grep("^sex$", names(poplong), value=TRUE, ignore.case=TRUE)


if("GROUP" %in% names(poplong)) {
    if(!length(sexcol)){
 # n <- regexpr("_", poplong$GROUP, fixed = TRUE)
 # poplong$AGE = substr(poplong$GROUP, n-1, 100)
#  poplong$AGE = substr(poplong$GROUP, 3, 4)
  poplong$sex = factor(substr(poplong$GROUP, 1, 1))}
     if(!length(agecol)){
   ageNumeric = as.numeric(substr(poplong$GROUP, 3, 4))
   poplong$age = cut(ageNumeric, ageBreaks$breaks, right=F)
  #Get rid of M/F if age group has only one digit
#  ageterm <- c(grep("^M", poplong$AGE, value=TRUE), grep("^F", poplong$AGE, value=TRUE))
#  agetermIndex <- c(grep("^M", poplong$AGE), grep("^F", poplong$AGE))
#  poplong$AGE[agetermIndex]<- substr(ageterm,2,100)
  }else {
  warning("no age and sex variables found or no group variable found in popdata")
  }
}


# aggregate if necessary
#if(!breaks %in% ageBreaks$breaks){

if(ageBreaks$mustAggregate == TRUE & is.null(aggregate.by) ){
popa = aggregate(poplong$POPULATION, list(age= poplong$age), sum)
names(popa)[names(popa)=="x"] = "POPULATION"
poplong <- popa
}



#attributes(poplong)$ageBreaks = ageBreaks$breaks
#return(poplong)



poplong$id<-NULL
row.names(poplong)<-NULL
#names(poplong)<-toupper(names(poplong))

if(!is.null(aggregate.by)) {

  popa <- aggregate(poplong$POPULATION, poplong[, aggregate.by, drop=F], sum)

  # change x column name to 'population'
  names(popa)[names(popa)=="x"] = "POPULATION"
  names(popa)[names(popa)=="poplong[, aggregate.by]"] = aggregate.by
  poplong<-popa

}

attributes(poplong)$breaks = ageBreaks

poplong


}
