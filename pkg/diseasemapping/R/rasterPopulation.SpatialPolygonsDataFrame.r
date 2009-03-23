
rasterPopulation.SpatialPolygonsDataFrame<-function(popdata,bbox=NULL,nrows=200, ncols=200, xmn=NULL, xmx=NULL, ymn=NULL, ymx=NULL, projs="NA"){


#if bbox is given, use it to create the raster,
if(!is.null(bbox)) {
  r <- raster(bbox,nrows=nrows,ncols=ncols,projs=projs)
}else if(all(!is.null(xmn),!is.null(xmx),!is.null(ymn),!is.null(ymx))){ #else if bbox is supplied manually
  r <- raster(nrows=nrows,ncols=ncols,xmn=xmn, xmx=xmx, ymn=ymn, ymx=ymx, projs=projs)
}else{ #else do not use a bbox, use sp's own bbox
  bbox<-extent(matrix(popdata@bbox,nrow=2))
  r <- raster(bbox,nrows=nrows,ncols=ncols,projs=projs)
}

#Grep the field we need
field<-grep("^expected",names(popdata),ignore.case = TRUE)
fieldname<-grep("^expected",names(popdata),ignore.case = TRUE,value=T)
r <- polygonsToRaster(popdata, r,field=field)

#convert to a sp dataframe, NA excluded
sp<-rasterToPoints(r,asSpatialPoints=T)
names(sp)<-fieldname
sp

}