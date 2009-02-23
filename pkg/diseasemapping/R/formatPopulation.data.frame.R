formatPopulation <- function(popdata, aggregate.by, breaks=NULL, years = NULL, ...) {
  UseMethod("formatPopulation")
}

`formatPopulation.data.frame` <-
function(popdata, aggregate.by=NULL, years=NULL, breaks=NULL, ...) {

#popdata <- popdata@data
ageBreaks = getBreaks(names(popdata), breaks)

####reshape the popdata:
poplong = reshape(popdata,  varying=ageBreaks$oldNames, direction="long",
	v.names="POPULATION", timevar="GROUP", times = ageBreaks$newNames)
# create age and sex variables
age = as.numeric(substr(poplong$GROUP, 3, 4))
poplong$cutAge = cut(age, ageBreaks$breaks, right=F)
sex = substr(poplong$GROUP, 1, 1)	

# aggregate if necessary
#if(!breaks %in% ageBreaks$breaks){
 if( ageBreaks$mustAggregate == TRUE & is.null(aggregate.by) ){
popa = aggregate(poplong$POPULATION, list(age=poplong$cutAge), sum)
names(popa)[names(popa)=="x"] = "POPULATION"
poplong <- popa
}

if(length(years)) {
	Ndata = dim(popdata)[1]
	popdata = popdata[rep(1:Ndata, rep(length(years), Ndata)), ]
	popdata$year = rep(years, Ndata)
}


#attributes(poplong)$ageBreaks = ageBreaks$breaks
#return(poplong)



if(! all(c("AGE", "SEX") %in% names(poplong))) {
  if("GROUP" %in% names(poplong)) {

 # n <- regexpr("_", poplong$GROUP, fixed = TRUE)
 # poplong$AGE = substr(poplong$GROUP, n-1, 100)

  poplong$AGE = substr(poplong$GROUP, 3, 4)
  poplong$SEX = factor(substr(poplong$GROUP, 1, 1))
  #Get rid of M/F if age group has only one digit
  ageterm <- c(grep("^M", poplong$AGE, value=TRUE), grep("^F", poplong$AGE, value=TRUE)) 
  agetermIndex <- c(grep("^M", poplong$AGE), grep("^F", poplong$AGE))
  poplong$AGE[agetermIndex]<- substr(ageterm,2,100)  
  }else {
  warning("no age and sex variables found or no group variable found in popdata")
  }
}

poplong$id<-NULL
row.names(poplong)<-NULL
names(poplong)<-toupper(names(poplong))

if(!is.null(aggregate.by)) {

  popa <- aggregate(poplong$POPULATION,poplong[,toupper(aggregate.by),drop=FALSE], sum)
  # change x column name to 'population'
  names(popa)[names(popa)=="x"] = "POPULATION"
  poplong<-popa

}

attributes(poplong)$breaks = ageBreaks

poplong


}
