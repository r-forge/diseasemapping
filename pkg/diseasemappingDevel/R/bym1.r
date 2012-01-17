`bym` = function(formula="observed ~ 1", family = "poisson",offset="logExpected", spPolygon=sp,ID="ID",
      prior.Spatial = c(1,0.001), prior.iid=c(1,0.001), debug = FALSE,verbose=TRUE,keep=TRUE,exceed=1.5){
      
        newdata = spPolygon@data
        nb = poly2nb(spPolygon)
	#create graph file
	nb2inla("graph.graph",nb)
 
 
      
        param = paste(c(prior.iid,prior.Spatial),collapse=",")
        formula.byn = as.formula(paste(formula,"+","f(",ID,",model=\"bym\",graph.file=\"graph.graph\",param=c("      				,param,"),values=",ID,")",sep=""))
      
       
        if(is.character(offset)) offset = data[,offset]
      

  
    result1.bym = inla(formula.bym,family=family, data=newdata, offset=offset, control.results=list(return.marginals.random=TRUE), control.predictor = list(compute = TRUE),verbose=verbose,debug=debug,keep=TRUE)

    cat("Colletcting results...\n")
    n<-dim(newdata)[1]
    newdata$BYM = result1.bym$summary.random[[1]][1:n,2]
    newdata$SPTIAL = result1.bym$summary.random[[1]][(n+1):(2*n),2]
    #newdata$fitted = result1.bym$summary.fitted.values[,1]/newdata$expcted
    #newdata$linearWoff = result1.bym$summary.linear.predictor[,1]
    newdata$linear =  result1.bym$summary.linear.predictor[,1] - newdata$logExpected

    newdata$fittedWoff= 0

    for(i in 1:dim(newdata)[1]){ ##by inla
    newdata$fittedWoff[i] = inla.expectation(exp, result1.bym$marginals.linear.predictor[[i]])}
    newdata$fitted = newdata$fittedWoff/newdata$expected

    #Exceedance
    MarginalBYM<-result1.bym$marginals.random[[1]][1:n]
    MarginalSPATIAL<-result1.bym$marginals.random[[1]][(n+1):(2*n)]

   
    newdata[,paste("Uex",1.5,sep="")]<-getEx(MarginalU,exceed)
    newdata[,paste("fitted",1.5,sep="")]<-getExOFF(result1.bym$marginals.linear.predictor,exceed,offsetExp = newdata$expected)
    


    spPolygon@data = newdata
    m = list(spPolygon, result1.bym) ;names(m) = c("Polygon","INLA")
    m
}


