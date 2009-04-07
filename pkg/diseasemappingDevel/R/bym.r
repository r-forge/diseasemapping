`bym` = function(formula="observed ~ 1", family = "poisson", data,offset="logExpected", spatial="region", graph.file="graph",
      control.Spatial = "param=c(1,0.01),initial=4", control.iid="param=c(1,0.001),initial=4", debug = FALSE,verbose=TRUE,keep=TRUE){
      
      
      if(length(grep("besag",formula,value=T)) == 0){
      
        spatial.name<-paste(spatial,"SP",sep="")
        data[,spatial.name]<-data[,spatial]
        formula1 = as.formula(paste(formula,"+","f(",spatial.name,",model=","\"besag\",","graph.file=\"",graph.file,"\",",control.Spatial,")",
                  "+","f(",spatial,",model=","\"iid\",",control.Spatial,")",sep=""))
      
      
      }



    inla(formula=formula1, family=family, data=data, offset=get(offset), control.results=list(return.marginals.random=TRUE), control.predictor = list(return.marginals=TRUE, compute = TRUE),
          verbose=verbose, midFunction=BYMmidFunction, debug=debug,keep=TRUE)
          

}


`BYMmidFunction` = function(inla.dir,n,struct,unstruct) {
    cat("Writting Linear Combination Files to", inla.dir,"\n")
    
     for (i in seq(0,n-1) ){
      inla.lincomb.section(inla.dir, lincomb.spec=16,num=i,struct=struct,unstruct=unstruct)
    }
    cat("Done...\n")
}

