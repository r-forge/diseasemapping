`inla.user.hook` = function(file.ini = NULL, data.dir = NULL, results.dir = NULL, formula = NULL, data = NULL)
{
    
    

     struct=strsplit(gsub("^[[:space:]]*f[[:space:]]*\\([[:space:]]*", "", grep("besag",strsplit(as.character(formula)[3],"\\+")[[1]],value=T)), ",")[[1]][1]
     
     unstruct=strsplit(gsub("^[[:space:]]*f[[:space:]]*\\([[:space:]]*", "", grep("iid",strsplit(as.character(formula)[3],"\\+")[[1]],value=T)), ",")[[1]][1]
    
        
     n<-dim(data)[1]
     cat(paste("Modifying",file.ini,sep=""),"\n", sep = " ")
     cat(paste("Writting Linear Combination data files to ",data.dir,sep=""),"\n", sep = " ")

     for (i in seq(0,n-1) ){
      inla.lincomb.section(file.ini, data.dir,lincomb.spec=16,num=i,struct=struct,unstruct=unstruct)
    }

    cat("Done...\n")


}


`inla.lincomb.section`=function (file.ini, data.dir,lincomb.spec,num,struct,unstruct)
{

     num<-inla.num(num)
    ############Modify ini file

    cat(paste("[Lincomb",num,"]\n",sep=""), sep = " ", file = file.ini, append = TRUE)
    cat("type = lincomb\n", sep = " ", file = file.ini, append = TRUE)
    cat(paste("dir","=","fixed.effect","9999",num,"\n",sep=""),file = file.ini, append = TRUE)
    cat("filename = ", paste("$DATADIR/lincomb",num,".dat",sep=""), "\n", sep = " ", file = file.ini,append = TRUE)
    if (!is.null(lincomb.spec)) {
        cat("precision = ", lincomb.spec, sep = " ", file = file.ini,
            append = TRUE)
        cat("\n", sep = " ", file = file.ini, append = TRUE)
    }
   cat("\n", sep = " ", file = file.ini, append = TRUE)

    ############Write Data files

    datafile<-paste(data.dir,paste("lincomb",num,".dat",sep=""),sep="/")

      cat(struct, sep = " ", file = datafile)
        cat("\n", sep = " ", file = datafile, append = TRUE)
      cat(paste(num,1.0,sep=" "), sep = " ", file = datafile,append = TRUE)
          cat("\n", sep = " ", file = datafile, append = TRUE)
          cat("\n", sep = " ", file = datafile, append = TRUE)
      cat(unstruct, sep = " ", file = datafile,append = TRUE)
          cat("\n", sep = " ", file = datafile, append = TRUE)
          cat(paste(num,1.0,sep=" "), sep = " ", file = datafile,append = TRUE)
      cat("\n", sep = " ", file = datafile, append = TRUE)

}
