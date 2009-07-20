getINLAData<-function(coarseRaster,casedata){

 casedata$one <- 1
 if(class(casedata)=="SpatiaPointsDataFrame") casedata$coords<-casedata@coords
 if(is.null(casedata$coords)) warning("No exact locations for case data")

 expected<-raster2DF(coarseRaster)
 names(expected)[3] = "expected"
 observed <- data.frame(rasterToPoints(pointsToRaster(coarseRaster,xy=casedata$coords,val=casedata$one)))
 names(observed)[3] = "observed"

 all<-merge(expected, observed,all.x=T)
 all$observed[is.na(all$observed)] = 0
 all$observed[is.na(all$expected)] = NA
 all$observed[all$expected==0] = NA
 all$E = all$expected*xres(coarseRaster)*yres(coarseRaster)
 all$logE = log(all$E)
 all$logE[all$logE == -Inf] = 0 #this could be any number
 all$logE[is.na(all$logE)] = 0

 all$ID<-1:dim(all)[1]

 all
}


raster2DF<-function(raster){
   data.frame(x=xyFromCell(raster,1:(attributes(raster)$nrow*attributes(raster)$ncols))[,1],
    y=xyFromCell(raster,1:(attributes(raster)$nrow*attributes(raster)$ncols))[,2],
    value = values(raster))
}
