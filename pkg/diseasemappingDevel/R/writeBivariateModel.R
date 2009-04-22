
writeBivariateModel <- function() { # need to write in the arguments!


  sink(file)
  
  cat("model{\n\n")

# write smoking model
writeBugsModel(prefix="", file=NULL)

# write cancer model

writeBugsModel(prefix="smoking", file=NULL)


# put smoking in cancer!


  cat("\n} # model\n") 

  sink()


}