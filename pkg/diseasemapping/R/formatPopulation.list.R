`formatPopulation.list`<-function(popdata, aggregate.by=NULL, cyears=as.integer(names(popdata)),  times="YEAR", ...) {
    
    times<-toupper(times)
    
   #If aggregate, then see if YEAR is there or not, if so, remove it
   if(!is.null(aggregate.by)){ 
   agg<-aggregate.by<-toupper(aggregate.by)
   byYear<- times %in% aggregate.by
   if(byYear){aggregate.by<-aggregate.by[-which(aggregate.by==times)]}
   }

   #if(class(popdata[[1]])== "SpatialPolygonsDataFrame"){
   #listpop<-lapply(popdata,formatPopulation.SpatialPolygonsDataFrame,aggregate.by)
   #}else{
   #listpop<-lapply(popdata,formatPopulation.data.frame,aggregate.by)
   #}
  
   listpop<-lapply(popdata,formatPopulation,aggregate.by)
  
   breaks = attributes(listpop[[1]])$breaks
  
   listdataframe<-lapply(listpop,as.data.frame)
   pop<-NULL
   for (i in 1:length(listdataframe)){
    temp<-listdataframe[[i]]
    temp[,times]<-cyears[i]
    pop<-rbind(pop,temp)
   }
   
   attributes(pop)$breaks = breaks
  
   
   pop
   #if(byYear) {p<- aggregate(pop$POPULATION,pop[,agg,drop=FALSE],sum)}
   #names(p)[names(p)=="x"] = "POPULATION"
   #p
}