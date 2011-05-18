#popdata, a sp object
#given a bounding box or coordinates, otherwise using bbox from popdata
#return a stack of rasters with values specified in columns argument
#if cellCoarse and fact is null, then fine raster will be returned, otherwise the Coarse raster with fact
#cellFine and cellCoarse are dimention of cells in the same unit as the unit of popdata

rasterSMR.SpatialPolygonsDataFrame<-function(popdata,bbox=NULL,
		cellFine=c(50,50), fact=NULL, cellCoarse=NULL,
  xmn=NULL, xmx=NULL, ymn=NULL, ymx=NULL, 
  columns=c("expected_sqk"),proj4string=popdata@proj4string){

#aggregation indicator, if either one is not NULL, aggregate
#agg<-(!is.null(fact) | !is.null(cellCoarse))
cellCoarse=cellFine*fact
#if bbox is given, use it to create the raster,
if(!is.null(bbox)) {
  xmn<-xmin(bbox); xmx<-xmax(bbox); ymn<-ymin(bbox); ymx<-ymax(bbox);
#else if coords are not supplies fully, use sp's own bbox
}else if(any(is.null(xmn),is.null(xmx),is.null(ymn),is.null(ymx))){
  bbox<-extent(matrix(popdata@bbox,nrow=2))
  xmn<-xmin(bbox); xmx<-xmax(bbox); ymn<-ymin(bbox); ymx<-ymax(bbox);
}

#if no aggregation needed, adjust bbox according to fine cells
if(length(cellCoarse)==0){cellCoarse<-cellFine}

#Make sure bbox is multiples of coarse, otherwise ajusted bbox
xdis<-xmx-xmn; ydis<-ymx-ymn
xmax<-xmx;xmin<-xmn;ymax<-ymx;ymin<-ymn
#get reminders
x=xdis %% cellCoarse[1]; y=ydis %% cellCoarse[2];
if(x!=0){xmax<-xmx + (cellCoarse[1]-x)/2; xmin<-xmn - (cellCoarse[1]-x)/2
        warning("Bounding Box is adjusted on X direction to make it multiples of coarse cells: ",paste(xmin,xmax,sep=","))}
if(y!=0){ymax<-ymx + (cellCoarse[2]-y)/2; ymin<-ymn - (cellCoarse[2]-y)/2
        warning("Bounding Box is adjusted on Y direction to make it multiples of coarse cells: ",paste(ymin,ymax,sep=","))}
        
        
#check if Coarse is mutilples of fine
#if not, adjust fine and warning
if(!is.null(cellCoarse) & is.null(fact)){
cx<-cellCoarse[1] %% cellFine[1]
cy<-cellCoarse[1] %% cellFine[1]

#find aggregation factor
fact<-c(cellCoarse[1] %/% cellFine[1],cellCoarse[2] %/% cellFine[2])

if(cx!=0){ cellFine[1] <- cellCoarse[1]/fact[1]; warning("Fine cell on X direction has been adjusted to ",cellFine[1])}
if(cy!=0){ cellFine[2] <- cellCoarse[2]/fact[2]; warning("Fine cell on Y direction has been adjusted to ",cellFine[2])}

}


##Find number of row and cols base on fine cells
ncols <-  (xmax-xmin) %/% cellFine[1];  nrows <-  (ymax-ymin) %/% cellFine[2]
r <- raster(nrows=nrows,ncols=ncols,xmn=xmin, xmx=xmax, 
		ymn=ymin, ymx=ymax, crs=proj4string)

#Grep the columns we need
for (i in 1:length(columns)){
  f<-which(names(popdata) == columns[i])
  k<-polygonsToRaster(popdata, r,field=f)
  #if need to be aggregated

  if(!all(cellFine==cellCoarse)){
    k<-aggregate(k,fact,na.rm=T)
   }
   #make them a stack
  if(i==1){sta<-stack(k)}else{sta<-addRasters(sta,k)}
}

sta

}