
rasterSMR.SpatialPolygonsDataFrame<-function(popdata,bbox=NULL,nrows=200, ncols=200,
  xmn=NULL, xmx=NULL, ymn=NULL, ymx=NULL, columns=c("expected"),projs="NA"){


#if bbox is given, use it to create the raster,
if(!is.null(bbox)) {
  r <- raster(bbox,nrows=nrows,ncols=ncols,projs=projs)
}else if(all(!is.null(xmn),!is.null(xmx),!is.null(ymn),!is.null(ymx))){ #else if bbox is supplied manually
  r <- raster(nrows=nrows,ncols=ncols,xmn=xmn, xmx=xmx, ymn=ymn, ymx=ymx, projs=projs)
}else{ #else do not use a bbox, use sp's own bbox
  bbox<-extent(matrix(popdata@bbox,nrow=2))
  r <- raster(bbox,nrows=nrows,ncols=ncols,projs=projs)
}

#Grep the columns we need
field<-NULL;lr<-list()
for (i in 1:length(columns)){
  f<-which(names(popdata) == columns[i])
  k<-polygonsToRaster(popdata, r,field=f)
  lr[[i]]<-k
  field<-c(field,f)
}
names(lr)<-names(popdata)[field]

lr

}