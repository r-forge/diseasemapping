writePrior <- function(x, file="priors.txt") {

SpsTrans = NULL

  
  sink(file)
  cat("transition", "parameter", "distribution", "mean", "sd", sep="\t" )
  cat("\n")
  for(Dtrans in names(x)) {
    for(Dparam in names(x[[Dtrans]])) {
      distribution =   attributes(x[[Dtrans]][[Dparam]])$distribution
      if(distribution == "psPrior" ) {
        SpsTrans = cbind(SpsTrans, c(Dtrans, Dparam))  
      } else {
        cat(Dtrans, Dparam, distribution,
          x[[Dtrans]][[Dparam]][c("mean","sd")],  sep="\t")
        cat("\n")
      }
    }
  }
  
  if(ncol(SpsTrans)) {
    # names ofthe hyperparameters
    parNames = names(x[[SpsTrans[1,1] ]][[SpsTrans[2,1] ]])
    cat("PS Priors\n")
    cat("transition","parameter", "distribution", parNames, sep ="\t")
    cat("\n")
    for(D in 1:ncol(SpsTrans) ) {
        cat(SpsTrans[1,D], SpsTrans[2,D], "psPrior",
          unlist(x[[Dtrans]][[Dparam]][parNames]), sep="\t")
       cat("\n")
    }
  }
  
  
  sink()

}

readPrior = function(file="priors.txt") {
  x = try(read.table(file, header=T, as.is=T), silent=T)
  if(class(x)=="try-error") {
    psPos = scan(file, sep="\n",what="a",quiet=T)
    psPos = grep("^PS Priors[[:space:]]*", psPos)
    if(length(psPos) != 1) 
      warning("prior file ", file,  " isnt in the right format")
    x = read.table(file, header=T, as.is=T, nrow=psPos-2)
    psPriorTable = read.table(file, header=T, as.is=T, skip=psPos)
  } else {
    psPriorTable = NULL
  }
  
  
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
  
  # ps priors 
  if(!is.null(psPriorTable)) {
    parNames = names(psPriorTable)
    parNames = parNames[!parNames %in% c("transition","parameter","distribution")]
    for(Dtrans in unique(psPriorTable$transition)) {
      x[[Dtrans]] = list()
      thisPs = psPriorTable[psPriorTable$transition == Dtrans,]
      rownames(thisPs) = thisPs[,"parameter"]
      for(Dparameter in  unique(thisPs$parameter)) {
        x[[Dtrans]][[Dparameter]] = as.list(thisPs[Dparameter,parNames])
    }
  }
  }
  x
  
  


}

