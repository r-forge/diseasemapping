getPopDis<-function(popdata,location,digits=0,threshold=NULL,moveYear=NULL){

  #Find Center of each region
  for (i in 1:length(popdata)){
    popdata[[i]]$x.center<-NULL; popdata[[i]]$y.center<-NULL;
    popdata[[i]]$x.center<-lapply(popdata[[i]]@polygons,function(x) x@Polygons[[1]]@labpt[1])
    popdata[[i]]$y.center<-lapply(popdata[[i]]@polygons,function(x) x@Polygons[[1]]@labpt[2])
  }



 #####and create distance
 loc1<-complex(real=location[1],imaginary=location[2])


  for (i in 1:length(popdata)){
    cgrid<-complex(real=popdata[[i]]$x.center,imaginary=popdata[[i]]$y.center)
    popdata[[i]]$D<-round(as.vector(outer(cgrid,loc1,getd)),digits)
    #change distances bigger than threshold to the threshold
    if(!is.null(threshold)) popdata[[i]]$D[popdata[[i]]$D > threshold] = threshold
    if(!is.null(moveYear))popdata[[i]]$move<-as.factor(as.integer(names(popdata)[i]) > moveYear)

 }
 popdata
}

getd<-function(p1,p2){
    sqrt((Re(p1)-Re(p2))^2 + (Im(p1)-Im(p2))^2)
}