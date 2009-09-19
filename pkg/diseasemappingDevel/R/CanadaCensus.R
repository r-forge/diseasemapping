CanadaCensusCSD = function(csd="3507008") {



startString = "http://www12.statcan.gc.ca/census-recensement/2006/dp-pd/prof/92-591/details/download-telecharger/CSV.cfm?Lang=E&Geo1=CSD&Code1="

endString="&Geo2=PR&Code2=35&Data=Count&SearchText=&SearchType=&SearchPR=&B1=Families%20and%20households&Custom="


downloadstring = paste(startString,  csd, endString, sep="")

fromWeb = read.table(downloadstring, header=T, sep=",", skip=1, as.is=T, nrow=58, fill=T)[,-1]

incCol = grep("median income[[:print:]]+all[[:print:]]+household", fromWeb[,1], ignore.case=T)

if(length(incCol)) {
return(fromWeb[incCol,2])
} else {
 return(NA)
}


}


addIncome = function(x, verbose=F) {

  if("CSDUID" %in% names(x)) {
    x$medIncome = NA
    for(D in seq(1, length(x$CSDUID)) ) {
      Dchar = as.character(x$CSDUID[D])
if(verbose){
cat(Dchar, " ")
if((D/10) == round(D/10)) cat("\n")
}
        x$medIncome[D] = 
          CanadaCensusCSD(csd=Dchar)
    }  
if(verbose){
cat("\n")
}

  } else {
    warning("no csduid")
  }

   x

}