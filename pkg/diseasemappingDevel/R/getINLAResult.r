getINLAresult<-function(result,all,exceed=c(1.1,1.2),below=NULL,write.shape="shapefile",CRSString="+proj=utm +zone=17 +datum=NAD83"){

 all$U<-result$summary.random[[1]][,2]
 all$fitted<-result$summary.fitted.values[,1]
 all$linear<-result$summary.linear.predictor[,1]

 for (i in 1:length(exceed)){
 	fitName <-paste("fitted",exceed[i],sep="ex")
 	UName <-paste("U",exceed[i],sep="ex")
 	all[[fitName]]<-getEx(result$marginals.linear.predictor,exceed[i])
 	all[[UName]]<-getEx(result$marginals.random[[1]],exceed[i])
 }

 if(!is.null(below)){
 	for (i in 1:length(below)){
 		fitName <-paste("fitted",below[i],sep="_")
 		UName <-paste("U",below[i],sep="_")
 		all[[fitName]]<-getExx(result$marginals.linear.predictor,below[i])
 		all[[UName]]<-getExx(result$marginals.random[[1]],below[i])
	 }
 }

 if(length(write.shape)==1){
        library(maptools)
	sppp<-SpatialPointsDataFrame(all[,c("x","y")],all, proj=CRS(CRSString))
	writePointsShape(sppp, write.shape)
 }


all  
}


getEx<-function(post,r){
        V<-NULL
	for(i in 1:length(post)){
		post1<-post[[i]]
		colnames(post1) = c("x", "y")
                thediff = diff(post1[,"x"])
                thediff = 0.5*c(0, thediff) + 0.5*c(thediff, 0)
                above= post1[,"x"] > log(r)
		V[i]<-sum(post1[above,"y"] * thediff[above])
	}
       V

}


getExx<-function(post,r){
        V<-NULL
	for(i in 1:length(post)){
		post1<-post[[i]]
		colnames(post1) = c("x", "y")
                thediff = diff(post1[,"x"])
                thediff = 0.5*c(0, thediff) + 0.5*c(thediff, 0)
                above= post1[,"x"] < log(r)
		V[i]<-sum(post1[above,"y"] * thediff[above])
	}
       V

}
