#addCensus indicates if adding the values specified in columns argument over census
#return a list of rasters stacks if addCensus is False, otherwise a single stack
##Add up two vectors by element
add<-function(x,y){
xna = is.na(x)
yna = is.na(y)
result = x + y

onlyY = yna & !xna
result[onlyY] = x[onlyY]

onlyX = !yna & xna
result[onlyX] = y[onlyX]

result
}



rasterSMR.list<-function(poplist,bbox=NULL,xmn=NULL, xmx=NULL, ymn=NULL, ymx=NULL,
        cellFine=c(50,50), cellCoarse=NULL, fact=NULL,projs="NA",columns=c("expected_sqk"),addCensus=TRUE){

#if coordinates are not given, use bbox
if(any(is.null(xmn),is.null(xmx),is.null(ymn),is.null(ymx))){
#find a bbox for all census
x<-NULL;y<-NULL

for (year in 1:length(poplist)){
  x<-c(x,poplist[[year]]@bbox[1,])
  y<-c(y,poplist[[year]]@bbox[2,])
}

xlimits = c(min(x),max(x))
ylimits = c(min(y),max(y))
bbox<-extent(matrix(c(xlimits,ylimits),nrow=2,byrow=T))
}


a<-lapply(poplist,rasterSMR,bbox=bbox,fact=fact,cellFine=cellFine, cellCoarse=cellCoarse,xmn=xmn, xmx=xmx, ymn=ymn, ymx=ymx, projs=projs,columns=columns)

if(addCensus){

    #add up columns in each census, create a single stack
    for(j in 1:length(columns)){
      news<-a[[1]]@layers[[j]] #initialize
      for (i in 2:length(a)){
        news<-overlay(news,a[[i]]@layers[[j]],fun=add)
      }
     if(j==1){sta<-stack(news)}else{sta<-addRasters(sta,news)}
    }

}
if(addCensus) list(indivdidual=a,total=sta)
else a

}

