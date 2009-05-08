writePrior <- function(x, file="priors.txt") {
  
  sink(file)
  cat("transition", "parameter", "mean", "sd", "distribution\n")
  for(Dtrans in names(x)) {
    for(Dparam in names(x[[Dtrans]]))
      cat(Dtrans, Dparam, 
        x[[Dtrans]][[Dparam]][c("mean","sd")], 
      attributes(x[[Dtrans]][[Dparam]])$distribution,
      "\n", sep="\t")
  }
  
  sink()

}

readPrior = function(file="priors.txt") {
  x = read.table(file, header=T, as.is=T)
  
  x=by(x, x$transition, function(qq) {
      thisTrans = by(qq, 
            qq$parameter, function(ww) {
         thisDist = unlist(
              ww[1,c("mean","sd")])
         attributes(thisDist)$distribution <- 
              ww[,"distribution"]
         thisDist
     } )
     class(thisTrans)="list"
     thenames = names(thisTrans)
     attributes(thisTrans) = NULL
     names(thisTrans) = thenames
     priorShapeScale(thisTrans)
  } )
  
  class(x) = "list"
  thenames = names(x)
  attributes(x) = NULL
  names(x) = thenames
  x
  


}