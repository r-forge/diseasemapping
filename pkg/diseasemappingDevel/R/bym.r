`bym` = function(formula="observed ~ 1", family = "poisson", data,offset="logExpected", spatial="region", graph.file="graph",
      control.Spatial = "param=c(1,0.01),initial=4", control.iid="param=c(1,0.001),initial=4", debug = FALSE,verbose=TRUE,keep=TRUE){
      
 
 
      
      if(length(grep("besag",formula,value=T)) == 0){
      
        spatial.name<-paste(spatial,"SP",sep="")
        data[,spatial.name]<-data[,spatial]
        formula1 = as.formula(paste(formula,"+","f(",spatial.name,",model=","\"besag\",","graph.file=\"",graph.file,"\",",control.Spatial,")",
                  "+","f(",spatial,",model=","\"iid\",",control.Spatial,")",sep=""))
      
      
      }
      print(formula1)

   theHook = function(file.ini = NULL, data.dir = NULL, results.dir = NULL, 
    formula = NULL, data = NULL,args = NULL) {
    
         inla.lincomb.section(file.ini, data.dir,lincomb.spec=10,num=0,
          struct=spatial.name,unstruct=spatial,D=NULL)

    }


   if(is.character(offset))  {
      offset = data[,offset]
      }
    inla(formula=formula1, family=family, data=data, offset=offset, 
     control.results=list(return.marginals.random=TRUE), 
     control.predictor = list(compute = TRUE), # return.marginals=TRUE, 
     verbose=verbose,  debug=debug,keep=TRUE,
     user.hook=theHook)
          

}



