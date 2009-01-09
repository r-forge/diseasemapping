formatPopulation <- function(popdata, aggregate.by,breaks=NULL, ...) {
  UseMethod("formatPopulation")
}

`formatPopulation.data.frame` <-
function(popdata, aggregate.by=NULL, breaks=NULL, ...) {

#popdata <- popdata@data
ageBreaks = getBreaks(names(popdata))

####reshape the popdata:
poplong = reshape(popdata,  varying=ageBreaks$oldNames, direction="long",
	v.names="POPULATION", timevar="GROUP", times = ageBreaks$newNames)
	
attributes(poplong)$breaks = ageBreaks$something	
	
# check to see if the provided breaks are good
######## this is the part have the problem. 

    if (length(breaks) > 0) {      
    currentbreaks <- attributes(poplong)$breaks    
         if (length(currentbreaks) == 0){ 
              currentbreaks <- breaks}
         if (all(breaks %in% currentbreaks)) {
        attributes(poplong)$breaks = breaks
        } else {
        stop("population age breaks aren't nested in the breaks provided")}
    } else {breaks <- attributes(poplong)$breaks}
    if (length(breaks) == 0){
    breaks <- ageBreaks
    attributes(poplong)$breaks <- breaks
    }

# I used loop here, but took a long time, there should be an easy way to solve the problem.
# apply cut function with the desired breaks

#sex=factor(substr(poplong$GROUP, 1, 1))
a <- numeric()
age <- substr(poplong$GROUP, 3, 4)
for(i in 1:length(age)){
a[i] <- which(age[i] == unique(ageBreaks$age))
poplong$cutAge[i] <- ageBreaks$names[a[i]]
}

# aggregate if necessary
if(!breaks %in% ageBreaks){
popa <- aggregate(poplong$POPULATION, list(age=poplong$cutAge), sum)
names(popa)[names(popa)=="x"] = "POPULATION"
poplong <- popa
}else{poplong <- poplong}

#attributes(poplong)$ageBreaks = ageBreaks$breaks
#return(poplong)

 
#### need to change M85plus to "M85_89" for the reshape step:
if( all( c(grep("^M.*PLUS", names(popdata), value=TRUE),
    grep("^F.*PLUS", names(popdata), value=TRUE))  %in%
  names(popdata)) ){
a1 <- substr(names(popdata)[grep("^M.*PLUS", names(popdata))], 2, 3)
b1 <- paste(a1, as.numeric(a1)+4, sep="_")
c1 <- paste("M", b1, sep="")
a2 <- substr(names(popdata)[grep("^F.*PLUS", names(popdata))], 2, 3)
b2 <- paste(a1, as.numeric(a2)+4, sep="_")
c2 <- paste("F", b2, sep="")
names(popdata)[grep("^M.*PLUS", names(popdata))] <- c1
names(popdata)[grep("^F.*PLUS", names(popdata))] <- c2
}
thevarying = c(grep("^M.*_", names(popdata), value=TRUE), grep("^F.*_", names(popdata), value=TRUE))       



if(! all(c("AGE", "SEX") %in% names(poplong))) {
  if("GROUP" %in% names(poplong)) {

  n <- regexpr("_", poplong$GROUP, fixed = TRUE)
  poplong$AGE = substr(poplong$GROUP, n-2, 100)
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

if(is.null(aggregate.by) & !is.null(breaks) ) {
  # do some aggregation

}

poplong
}
