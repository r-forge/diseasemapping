`inla.user.hook` = function(file.ini = NULL, data.dir = NULL, results.dir = NULL, formula = NULL, data = NULL)
{
    
    

     struct=strsplit(gsub("^[[:space:]]*f[[:space:]]*\\([[:space:]]*", "", grep("besag",strsplit(as.character(formula)[3],"\\+")[[1]],value=T)), ",")[[1]][1]
     
     unstruct=strsplit(gsub("^[[:space:]]*f[[:space:]]*\\([[:space:]]*", "", grep("iid",strsplit(as.character(formula)[3],"\\+")[[1]],value=T)), ",")[[1]][1]
    
        
     n<-dim(data)[1]
     cat(paste("Modifying",file.ini,sep=""),"\n", sep = " ")
     cat(paste("Writting Linear Combination data files to",data.dir,sep=""),"\n", sep = " ")

     for (i in seq(0,n-1) ){
      inla.lincomb.section(file.ini, data.dir,lincomb.spec=16,num=i,struct=struct,unstruct=unstruct)
    }

    cat("Done...\n")


}
