#' Script for running inla remote

#+ library
library("INLA")
#'

#+ iniFileLinux
if(.Platform$OS.type == 'unix')
  cat(
      'RemoteINLA="/home/patrick/R/x86_64-pc-linux-gnu-library/3.4/INLA/bin/linux/64bit/inla"',
      'RemoteUserHost="patrick@darjeeling.pbrown.ca"',
      'SSHCMD="ssh"',    
      'SCPCMD="scp"',
      sep='\n', 
      file="~/.inlarc")
#'

#+ iniFileWindows
if(.Platform$OS.type == 'windows') {
  INLA::inla.setOption("cygwin",
      "C:\\\\Users\\\\brownpa\\\\programs\\\\cygwin"
  )
  cat(    
      'RemoteINLA="/home/patrick/R/x86_64-pc-linux-gnu-library/3.4/INLA/bin/linux/64bit/inla"',
      'RemoteUserHost="patrick@darjeeling.pbrown.ca"',
      'SSHCMD="/cygdrive/c/Users/brownpa/programs/x2go/plink"',
      'SCPCMD="/cygdrive/c/Users/brownpa/programs/x2go/pscp"',
      sep='\n', 
      file=normalizePath(file.path(
              Sys.getenv("USERPROFILE"), 
              ".inlarc"), '/', FALSE))
}
#'


#+ copyBatchFile
inlaOrig = system.file(paste('bin/remote/inla.remote', c('', '.cygwin'), sep=''), package='INLA')
inlaOrigBak = paste(inlaOrig, '.bak', sep='')
for(D in 1:length(inlaOrig)) {
  file.copy(inlaOrig[D], inlaOrigBak[D], overwrite=FALSE)
  file.copy(
      system.file('src/inla.remote.simple', package='Pmisc'),
      inlaOrig[D],
      overwrite=TRUE
  )
}
#'

#+ runInla
res = inla(y ~ 1 + x, 
    data = data.frame(y=1:3,x=factor(c("a","b",NA))), 
    verbose=TRUE,
    inla.call="remote")
#'

#+ results
knitr::kable(res$summary.fixed, digits=4)
#'

#+ restoreBatchFile
for(D in 1:length(inlaOrig)) {
  file.copy(inlaOrigBak[D], inlaOrig[D], overwrite=TRUE)
}
#'
