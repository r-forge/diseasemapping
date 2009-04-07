`inla.lincomb.section`=function (dir, lincomb.spec,num,struct,unstruct)
{

     file = paste(dir,"Model.ini",sep="/")
    ############Modify ini file
    
    cat(paste("[Lincomb",num,"]\n",sep=""), sep = " ", file = file, append = TRUE)
    cat("type = lincomb\n", sep = " ", file = file, append = TRUE)
    cat(paste("dir","=","fixed.effect","9999",num,"\n",sep=""),file = file, append = TRUE)
    cat("filename = ", paste("$DATADIR/lincomb",num,".dat",sep=""), "\n", sep = " ", file = file,append = TRUE)
    if (!is.null(lincomb.spec)) {
        cat("precision = ", lincomb.spec, sep = " ", file = file,
            append = TRUE)
        cat("\n", sep = " ", file = file, append = TRUE)
    }
   cat("\n", sep = " ", file = file, append = TRUE)

    ############Write Data files

    ddir<-paste(dir, "data.files/",sep="/")
    datafile<-paste(ddir,paste("lincomb",num,".dat",sep=""),sep="/")
    
     

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
