#return a list of dataframes if multiple census are given
#Need case file has region in it

`getSMR` <-
function(model, popdata, casedata=NULL,regionCode="CSDUID",
  regionCodeCases="CSD2006",cyears=NULL,year.range=NULL,area=FALSE,
  area.scale=1) {

morethanoneyear = class(popdata)=="list"

#see if it is a SP or a dataframe
if(morethanoneyear){
  isSP = (class(popdata[[1]])== "SpatialPolygonsDataFrame")
}else{
  isSP = (class(popdata)== "SpatialPolygonsDataFrame")
}

#if years not supplied, use the names of list
if(is.null(cyears) & morethanoneyear ){
  cyears = as.integer(names(popdata))
} 


#if more than one year, reshape bind them into dataframe,esle reshape only
#if(morethanoneyear) {
#poplong<-bindPopulation(popdata,cyears) 
#}else{
#poplong <- formatPopulation(popdata)
#}


if(area & isSP){
  if(morethanoneyear){
     #a list of vectors
     areas<-lapply(lapply(popdata,area),as.numeric)
     for (i in 1:length(popdata)){popdata[[i]]$sqk<-areas[[i]]*area.scale}          
  }else{
     popdata$sqk<-c(unlist(area(popdata)))*area.scale
  } 
}



#if more than one year, reshape bind them into dataframe,esle reshape only, no aggregation
#if(morethanoneyear) {
#poplong<-formatPopulation.list(popdata,cyears) 
#}else{
 #call diff functions if SP
# if(isSP) {poplong <- formatPopulation.SpatialPolygonsDataFrame(popdata)
# }else{poplong <- formatPopulation.data.frame(popdata)}
#}

poplong <- formatPopulation(popdata)

if(is.null(year.range) & morethanoneyear){
  year.range = range(poplong$YEAR)
} 

# get rid of age and sex groups which are not in the model
#  -- because if they're not in the model they have zero rates
if(length(model$sexSubset)==1) {
  # the 'getRates' function puts a sexSubset element showing which sexes were
  #  included in the model
  poplong=poplong[poplong$sex==model$sexSubset,]
}

if (morethanoneyear){
####find Popoluation census year
#a vector of all years
times<-c(year.range[1],cyears,year.range[2])
inter<-diff(times)/2 #mid points
#sum of consective mid points
nseq<-1:length(inter)-1
mseq<-2:length(inter)
interval<-inter[mseq] + inter[nseq]

names(interval)<-names(table(poplong$YEAR))
poplong$yearsForCensus = interval[as.character(poplong$YEAR)]
poplong$POPULATION = poplong$POPULATION  * poplong$yearsForCensus 
poplong$YEAR= factor(poplong$YEAR, levels = unique(poplong$YEAR))  
poplong$YEAR = factor(poplong$YEAR)
}

poplong <- poplong [poplong$POPULATION>0,]
poplong$LOGPOP = log(poplong$POPULATION)

poplong$SEX = factor(poplong$SEX)
poplong$AGE = factor(poplong$AGE)

names(poplong)<-tolower(names(poplong))

# find the ages used in the model and get rid of those that do not exist in model
#modelAges = model$xlevels$age
#poplong = poplong[poplong$age %in% modelAges,]

#remove main effects that not in the model
for(Dlevel in names(model$xlevels)) {
alllevels = levels(poplong[[Dlevel]])
if(!all(alllevels %in% model$xlevels[[Dlevel]])) {
	# remove the rows from the missing levels
	tokeep = poplong[[Dlevel]] %in% model$xlevels[[Dlevel]]
	poplong = poplong[tokeep,]
 }
}


#Take out interactions with NA, work with model without interactions
#temp<-row.names(summary(model)$coefficients)[-1]
interactNA<-names(model$coefficients)[is.na(model$coefficients)]
#if there is any coeff with NA, remove
if(length(interactNA)>0){
  interact<-grep(":",interactNA,value=TRUE)
  #keep the rows that the interaction is not NA,
  poplong$param<-paste(paste("age",poplong$age,sep=""),paste("sex",poplong$sex,sep=""),sep=":")
  poplong = poplong[!poplong$param %in% interact,]
  poplong$param<-NULL
}
#find expected
if(morethanoneyear) {   
agg<-c("year","age","sex","logpop")
}else{
agg<-c("age","sex","logpop")
}
poplong$expected<-predict(model, poplong[,agg], type="response")

#fing expected per sqK
if(area){
poplong$expected_sqk<-poplong$expected/poplong$sqk
}

#aggregate by region 
regionCode<-tolower(regionCode)
if(morethanoneyear) {
  if(area){
  newpop<-aggregate(poplong[,c("expected","expected_sqk")],list(Region=poplong[[regionCode]],Year=poplong$year),sum)
  }else{newpop<-aggregate(poplong[,"expected"],list(Region=poplong[[regionCode]],Year=poplong$year),sum)
        names(newpop)[names(newpop)=="x"] = "expected"}
}else{
  if(area){
    newpop<-aggregate(poplong[,c("expected","expected_sqk")],
      list(Region=poplong[[regionCode]]),sum)
  }else{
    newpop<-aggregate(poplong[,"expected"],
      list(Region=poplong[[regionCode]]),sum)
      names(newpop)[names(newpop)=="x"] = "expected"
  }
}

#names(newpop)[names(newpop)=="x"] = "expected"
newpop$logExpected = log(newpop$expected)
#return this if case file is null
populationData<-newpop

############################################Case File
if(!is.null(casedata)) {
  # get the total number of cases per region per year, using the cases supplied
  casedata <- formatCases(casedata)
  # only use cases which are in the regions of the dataset
  casedata = casedata[as.character(casedata[,regionCodeCases]) %in%
    as.character(populationData[, "Region"]),]
  
  casecol = grep("^cases$", names(casedata), value=TRUE, ignore.case=TRUE)
  if(!length(casecol)) {
    #there is no case col
    casecol = "cases"
    cases[,casecol] = 1
  }

  
  if(morethanoneyear) {
  newcase<-aggregate(casedata[[casecol]],
    list(Region=casedata[[toupper(regionCodeCases)]],Year=casedata$YEAR),sum)
  }else{
  newcase<-aggregate(casedata[[casecol]],
    list(Region=casedata[[toupper(regionCodeCases)]]),sum)
  }
  names(newcase)[names(newcase)=="x"] = "Observed"

  #merge with population get full dataset, named population data
  populationData = merge(newpop, newcase,all.x=TRUE)
  populationData$Observed[is.na(populationData$Observed)] = 0
    
  
   if(!is.null(populationData$Observed)) {
  # if there is a 'cases' column, possibly computed from a case file (above)
  # then compute the SMR
  populationData$SMR <- populationData$Observed/populationData$expected
  }     

}#end of big if         


 #if more than one year and not SP object, create a list of df  
 if(morethanoneyear) {
   listpop<-list()
   for (i in 1:length(cyears)){
     listpop[[i]]<- merge(popdata[[i]],
      populationData[populationData$Year==cyears[i],],
      by.x=grep(regionCode, names(popdata[[i]]), ignore.case=TRUE, value=TRUE),
      by.y="Region", all.x=TRUE
     )
   }
   if(isSP) {
    for (i in 1:length(cyears)){
      popdata[[i]]@data = listpop[[i]]
    }
    listpop = popdata
   }
 } else {  # only one year
   listpop = merge(popdata@data, populationData,
      by.x=grep(paste("^", regionCode, "$", sep=""),
        names(popdata@data), ignore.case=TRUE, value=TRUE),
      by.y="Region", all.x=TRUE)
   if(isSP) {
    popdata@data = listpop
    listpop = popdata
   }
}
 
 
listpop
}


