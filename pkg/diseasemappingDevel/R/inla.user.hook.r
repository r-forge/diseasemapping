`inla.user.hook` = function(file.ini = NULL, data.dir = NULL, results.dir = NULL, formula = NULL, data = NULL,args = NULL)
{



     #struct=strsplit(gsub("^[[:space:]]*f[[:space:]]*\\([[:space:]]*", "", grep("besag",strsplit(as.character(formula)[3],"\\+")[[1]],value=T)), ",")[[1]][1]

     #unstruct=strsplit(gsub("^[[:space:]]*f[[:space:]]*\\([[:space:]]*", "", grep("iid",strsplit(as.character(formula)[3],"\\+")[[1]],value=T)), ",")[[1]][1]



     #n<-544
     cat(paste("Modifying ",file.ini,sep=""),"\n", sep = " ")
     cat(paste("Writting Linear Combination data files to ",data.dir,sep=""),"\n", sep = " ")


     #for (i in seq(0,n-1) ){

     inla.lincomb.section(file.ini, data.dir,lincomb.spec=10,num=0,struct="disCov",unstruct="disCov1",D="disCov2")

#200
#      inla.lincomb.section(file.ini, data.dir,lincomb.spec=10,num=1,struct="disCov",unstruct="disCov1",D="disCov2",wt1=0.85,wt2=0.96,wt3=0.99)


#500
#      inla.lincomb.section(file.ini, #data.dir,lincomb.spec=10,num=2,struct="disCov",unstruct="disCov1",D="disCov2",wt1=0.368,wt2=0.779,wt3=0.939)

#800
#      inla.lincomb.section(file.ini, #data.dir,lincomb.spec=10,num=3,struct="disCov",unstruct="disCov1",D="disCov2",wt1=0.077,wt2=0.53,wt3=0.852)

#1000
#      inla.lincomb.section(file.ini, #data.dir,lincomb.spec=10,num=4,struct="disCov",unstruct="disCov1",D="disCov2",wt1=0.018,wt2=0.368,wt3=0.779)

#1500
#      inla.lincomb.section(file.ini, #data.dir,lincomb.spec=10,num=5,struct="disCov",unstruct="disCov1",D="disCov2",wt1=0.0001,wt2=0.11,wt3=0.57)

#2000

#      inla.lincomb.section(file.ini, #data.dir,lincomb.spec=10,num=6,struct="disCov",unstruct="disCov1",D="disCov2",wt1=0.0000001,wt2=0.018,wt3=0.368)

#2500

 #     inla.lincomb.section(file.ini, data.dir,lincomb.spec=10,num=7,struct="disCov",unstruct="disCov1",D="disCov2",wt1=0,wt2=0.0019,wt3=0.21)

#3000

#      inla.lincomb.section(file.ini, data.dir,lincomb.spec=10,num=8,struct="disCov",unstruct="disCov1",D="disCov2",wt1=0,wt2=0.0001,wt3=0.11)
   # }

    cat("Done...\n")


}


`inla.lincomb.section`=function (file.ini, data.dir,lincomb.spec,num,struct,unstruct,D,wt1=1.0,wt2=1.0,wt3=1.0,index=0)
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
          cat(paste(index,wt1,sep=" "), sep = " ", file = datafile,append = TRUE)
          cat("\n", sep = " ", file = datafile, append = TRUE)
          cat("\n", sep = " ", file = datafile, append = TRUE)
      cat(unstruct, sep = " ", file = datafile,append = TRUE)
          cat("\n", sep = " ", file = datafile, append = TRUE)
          cat(paste(index,wt2,sep=" "), sep = " ", file = datafile,append = TRUE)
          cat("\n", sep = " ", file = datafile, append = TRUE)
          cat("\n", sep = " ", file = datafile, append = TRUE)
      cat(D, sep = " ", file = datafile,append = TRUE)
          cat("\n", sep = " ", file = datafile, append = TRUE)
          cat(paste(index,wt3,sep=" "), sep = " ", file = datafile,append = TRUE)
          cat("\n", sep = " ", file = datafile, append = TRUE)
          cat("\n", sep = " ", file = datafile, append = TRUE)

}
