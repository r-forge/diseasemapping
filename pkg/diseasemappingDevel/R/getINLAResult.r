getINLAResult<-function(result,all,exceed=c(1.1,1.2),below=NULL,write.shape="shapefile",CRSString="+proj=utm +zone=17 +datum=NAD83"){


 all$U<-result$summary.random[[1]][,2]
 all$fitted<-result$summary.fitted.values[,1]
    all$fitted[all$fitted<0]<-0
    all$linear<-result$summary.linear.predictor[,1]
    all$sdFitted = result$summary.fitted.values[,2]
    all$sdlinear<-result$summary.linear.predictor[,2]
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

 
 

 library(maptools)
 sppp<-SpatialPixelsDataFrame(all[,c("x","y")],all, proj=CRS(CRSString))
 
    if(length(write.shape)==1){
     	   writePointsShape(sppp, write.shape)
    }


sppp
}







getINLAResultCOV<-function(result,allPoly,partPoly=NULL,offsetExp=1,exceed=c(1.1,1.2),below=NULL){


 allPoly$U<-result$summary.random[[1]][,2]
 #all = allPoly@data

all = allPoly
 if(is.null(partPoly)){

    all$fitted<-result$summary.fitted.values[,1]/offsetExp
    all$fitted[all$fitted<0]<-0
    all$linear<-result$summary.linear.predictor[,1]-log(offsetExp)
    all$sdFitted = result$summary.fitted.values[,2]
    all$sdlinear<-result$summary.linear.predictor[,2]
      for (i in 1:length(exceed)){
 	      fitName <-paste("fitted",exceed[i],sep="ex")
 	      UName <-paste("U",exceed[i],sep="ex")
 	      all[[fitName]]<-getExOFF(result$marginals.linear.predictor,exceed[i],offsetExp)
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

 }


  if(!is.null(partPoly)){

    #part = partPoly@data
part = partPoly
    part$fitted<-result$summary.fitted.values[,1]/offsetExp
    part$fitted[part$fitted<0]<-0
    part$linear<-result$summary.linear.predictor[,1] - log(offsetExp)
    part$sdFitted = result$summary.fitted.values[,2]
    part$sdlinear<-result$summary.linear.predictor[,2]
      for (i in 1:length(exceed)){
 	      fitName <-paste("fitted",exceed[i],sep="ex")
        UName <-paste("U",exceed[i],sep="ex")
 	      part[[fitName]]<-getExOFF(result$marginals.linear.predictor,exceed[i],offsetExp)
 	      all[[UName]]<-getExOFF(result$marginals.random[[1]],exceed[i],offsetExp=NULL)
      }

    if(!is.null(below)){
 	    for (i in 1:length(below)){
 	    	fitName <-paste("fitted",below[i],sep="_")
 		    UName <-paste("U",below[i],sep="_")
 		    part[[fitName]]<-getExx(result$marginals.linear.predictor,below[i])
 		    all[[UName]]<-getExx(result$marginals.random[[1]],below[i])
	   }
    }
    #partPoly@data = part
   partPoly = part
 }

   #allPoly@data=all
allPoly=all

ll=list(allPoly,partPoly)
names(ll) = c("Grid","Partition")
ll
}




getExOFF<-function(post,r,offsetExp=NULL){
if(is.null(offsetExp)) offsetExp =rep(1,length(post))
if(length(offsetExp)==1 &offsetExp[1]==1) offsetExp =rep(1,length(post))
  offset = log(offsetExp)
  V<-NULL
	for(i in 1:length(post)){
		post1<-post[[i]]
		post1[,"x"] = post1[,"x"] - offset[i]
		colnames(post1) = c("x", "y")
                thediff = diff(post1[,"x"])
                thediff = 0.5*c(0, thediff) + 0.5*c(thediff, 0)
                above= post1[,"x"] > log(r)
		V[i]<-sum(post1[above,"y"] * thediff[above])
	}
       V

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
