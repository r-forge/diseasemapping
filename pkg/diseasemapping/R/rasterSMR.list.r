rasterSMR.list<-function(poplist,bbox=NULL,nrows=200, ncols=200, xmn=NULL, xmx=NULL, ymn=NULL, ymx=NULL, projs="NA",columns=c("expected")){


#find a bbox for all census
x<-NULL;y<-NULL

for (year in 1:length(poplist)){
  x<-c(x,poplist[[year]]@bbox[1,])
  y<-c(y,poplist[[year]]@bbox[2,])
}

xlimits = c(min(x),max(x))
ylimits = c(min(y),max(y))
bbox<-extent(matrix(c(xlimits,ylimits),nrow=2,byrow=T))

lapply(poplist,rasterSMR,bbox=bbox,nrows=nrows, ncols=ncols, xmn=xmn, xmx=xmx, ymn=ymn, ymx=ymx, projs=projs,columns=columns)


}