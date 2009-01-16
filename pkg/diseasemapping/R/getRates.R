#cyears: a vector of census years, default will use the names of popdata list
#year.range: A two-element numerical vector indicting the starting and ending year of study period, default set to range of casedata
#case.years: the variable name that stores the year variable of case file

`getRates` <-
function(casedata, popdata, formula, family=poisson, minimumAge=0,
   maximumAge=100, S=c("M", "F"),cyears=NULL,year.range=NULL,
   case.years=grep("^year$", names(casedata), ignore.case=TRUE, value=TRUE)[1], breaks=NULL){
# check the formula is one sided

      thevarying = c(grep("^M.*_", names(popdata), value=TRUE), grep("^F.*_", names(popdata), value=TRUE))       
      poplong = reshape(popdata,  varying=thevarying, direction="long", v.names="POPULATION", timevar="GROUP", times = thevarying)


    if (length(breaks) > 0) {      
    currentbreaks <- attributes(poplong)$breaks    
         if (length(currentbreaks) == 0){ 
              currentbreaks <- getBreaks(popdata)}
         if (all(breaks %in% currentbreaks)) {
        attributes(poplong)$breaks = breaks
        } else {
        stop("population age breaks aren't nested in the breaks provided")}
    } else {breaks <- attributes(poplong)$breaks}
    if (length(breaks) == 0){
    breaks <- getBreaks(popdata)
    attributes(poplong)$breaks <- breaks
    }


if(attributes(terms(formula))$response)
  warning("formula should be one sided")


morethanoneyear = class(popdata)=="list"

#is SP or not
if(morethanoneyear){
  isSP = (class(popdata[[1]])== "SpatialPolygonsDataFrame")
}else{
  isSP = (class(popdata)== "SpatialPolygonsDataFrame")
}

#if years not supplied, use the names of list
if(is.null(cyears) & morethanoneyear ){
  cyears = as.integer(names(popdata))
} 

#factors we need to aggregate by
theterms = (rownames(attributes(terms(formula))$factors))

#if more than one year, reshape bind them into dataframe,else reshape, the aggregate
#if(morethanoneyear) {
#pops<-formatPopulation.list(popdata,cyears,aggregate.by=theterms) 
#}else{
# #call diff functions if SP
# if(isSP) {pops <- formatPopulation.SpatialPolygonsDataFrame(popdata,aggregate.by=theterms)
# }else{pops <- formatPopulation.data.frame(popdata,aggregate.by=theterms)}
#}

pops <- formatPopulation(popdata,aggregate.by=theterms)

#format case data
casedata = formatCases(casedata)

#if ranges not supplied, use the year ranges of case files
if(is.null(year.range) & morethanoneyear){
  year.range = range(casedata[[case.years]])
} 



# keep only the desired sex in the dataset
if(length(S)==1) {
  if(length(grep("^sex$", theterms, ignore.case=TRUE)))
    warning("sex is in the model but only one sex is being used")
  casedata=casedata[casedata[[
    grep("^sex$", names(casedata), value=TRUE, ignore.case=TRUE)
      ]]==S,]
 }


#find number of cases per group
casecol = grep("^cases$", names(casedata), value=TRUE, ignore.case=TRUE)
if(!length(casecol)) {
  #there is no case col
  casecol = "cases"
  cases[,casecol] = 1
}

cases<-aggregate(casedata[[casecol]], casedata[,theterms], sum, na.rm=TRUE)
names(cases)[names(cases)=="x"] = "CASES"

#find population per group
#pops <- aggregate(poplong$POPULATION, poplong[,theterms,drop=FALSE], sum)
#names(pops)[names(pops)=="x"] = "POPULATION"

##### merge case data set and shape data set according to the same Year and DA2001
by.x =  paste("^", theterms, "$", sep="")
by.x = paste(by.x, collapse="|")
by.x = paste("(", by.x, ")", sep="")
by.pop = grep(by.x, names(pops), ignore.case=TRUE, value=TRUE)


newdata <- merge(cases, pops, by.x = theterms, by.y = by.pop)

if (morethanoneyear){
####find Popoluation census year
#a vector of all years
times<-c(year.range[1],cyears,year.range[2])
inter<-diff(times)/2 #mid points
#sum of consective mid points
nseq<-1:length(inter)-1
mseq<-2:length(inter)
interval<-inter[mseq] + inter[nseq]

names(interval)<-names(table(newdata$YEAR))
newdata$yearsForCensus = interval[as.character(newdata$YEAR)]
newdata$POPULATION = newdata$POPULATION  * newdata$yearsForCensus 
newdata$YEAR= factor(newdata$YEAR, levels = unique(newdata$YEAR))
}


newdata = newdata[newdata$POPULATION>0,]
newdata$logpop = log(newdata$POPULATION)
	

# make the age group with the most cases as the base line
agevar =  grep("^age$", theterms, ignore.case=TRUE)
if(length(agevar)==1) {
  agetable = tapply(newdata$CASES, newdata[[agevar]], sum)
  agetable = names(sort(agetable, decreasing=TRUE))
  newdata[[agevar]] = factor(as.character(newdata[[agevar]]),
    levels= agetable)
}


# add cases and logpop to formula
formula1 = update.formula(formula, CASES ~ offset(logpop) + .)
#return(newdata, formula1)

model = glm(formula1, family=family, data=newdata)

model$sexSubset = S

# model$years
# model$breaks

return(model)
}







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