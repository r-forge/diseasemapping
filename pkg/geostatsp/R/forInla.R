lcOneRow = function(thisrow, idxCol=NULL) {
  thisrow = thisrow[!is.na(thisrow)]
  if(length(thisrow)) {
    thisrow = sapply(thisrow, function(qq) list(list(weight=qq)))
    for(D  in idxCol)
      thisrow[[D]] = list(
        weight=1, 
        idx=thisrow[[D]]$weight
        )
    for(D in names(thisrow))
      thisrow[[D]] = thisrow[D]
    names(thisrow) = paste("v", 1:length(thisrow), sep="")
  }
  thisrow
}

inla.models=function(){
  if(requireNamespace("INLA", quietly=TRUE)){
    return(INLA::inla.models())
  } else {
    return(NULL)
  }
}

# TRUE only when the INLA package loads and its binaries work.
# requireNamespace alone is not enough: some installs have the R package
# but a missing/broken inla program, and inla.setOption then errors.
inlaAvailable = function() {
  if(!isTRUE(requireNamespace("INLA", quietly = TRUE)))
    return(FALSE)
  ok = try(INLA::inla.getOption("num.threads"), silent = TRUE)
  !inherits(ok, "try-error")
}

# Cap INLA threads when available; no-op (and no error) otherwise.
inlaSetThreads = function(num.threads = 2L) {
  if(!inlaAvailable())
    return(invisible(FALSE))
  try(INLA::inla.setOption(num.threads = num.threads), silent = TRUE)
  # not all versions of INLA support blas.num.threads
  try(INLA::inla.setOption(blas.num.threads = num.threads), silent = TRUE)
  invisible(TRUE)
}