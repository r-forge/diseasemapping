`inla.lincomb.section`=function (file.ini, data.dir,lincomb.spec,num,struct,unstruct)
{


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


