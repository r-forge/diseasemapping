formatPopulation <- function(popdata, aggregate.by, ...) {
  UseMethod("formatPopulation")
}

`formatPopulation.data.frame` <-
function(popdata,aggregate.by=NULL, ...) {

 
#if(class(popdata)== "SpatialPolygonsDataFrame"){
#pdataframe = popdata@data
#names(pdataframe)<-toupper(names(pdataframe))
#popdata <- pdataframe
#}

names(popdata)<-toupper(names(popdata))      

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


####reshape the popdata:
thevarying = c(grep("^M.*_", names(popdata), value=TRUE), grep("^F.*_", names(popdata), value=TRUE))       
poplong = reshape(popdata,  varying=thevarying, direction="long",
	v.names="POPULATION", timevar="GROUP", times = thevarying)

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
