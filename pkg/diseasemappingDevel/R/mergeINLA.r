mergeINLA<-function(data,INLA,spatial="regionSP",iid="region",exceed=1.2){

  data$spatial<-INLA$summary.random[[spatial]][,2]
  data$iid<-INLA$summary.random[[iid]][,2]
  data$fitted<-INLA$summary.fitted.values[,1]
  data$linear<-INLA$summary.fitted.values[,1]

  #Get Joint
  #temp<-data.frame(bym=INLA$summary.fixed[-1,1],ID=as.integer(substr(names(temp),8,100))+1)
  #t<-temp[order(temp$ID),]
  #data$bym<-t$bym
  data$bym<-INLA$summary.fixed[-1,1]


  #use predictos on the log scale to find out exceeding
  data$fitted_exceed<-getEx(INLA$marginals.linear.predictor,exceed)
  data$bym_exceed<-getEx(INLA$marginals.fixed[-1],exceed)


  data
}


#Computer exceedance proabilities

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
