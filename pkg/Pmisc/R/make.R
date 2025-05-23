#' @export
make = function(x, 
  suffix=NULL, beamer=FALSE,
  run=FALSE, ...) {

  if(!is.null(suffix))
    x = paste(
      tools::file_path_sans_ext(basename(x)),
      suffix, sep='.')
  
  theString= paste(' -f', 
    system.file('extsrc','knitrMakefile', package='Pmisc'),
    x)
  
  if(beamer)
    theString = paste(theString, 'pandocTo=beamer')
  
  dots = list(...)
  for(D in seq(from=1, by=1, len=length(dots)))
    theString = paste(theString, D)
  
  if(run) {
    system(paste('make ', theString))
  } else {
    cat(theString)
  }
  
  invisible(theString)
  
}

#' @export
Makefile = function(x, suffix,
  beamer = FALSE, output='Makefile', 
  cygwin=FALSE, fullpath=FALSE) {

  if(!missing(suffix) & !missing(x))
    x = paste(
      tools::file_path_sans_ext(basename(x)),
      suffix, sep='.')
  
  if(cygwin){
    whichfun = function(xx) {
      res = system(paste(
        Sys.which("bash"),
        " -c 'which --all ", xx, "'",
        sep=''
        ), intern=TRUE)
      if(length(res > 1)) {
        notCygdrive = grep("cygdrive", res, value=TRUE, invert=TRUE)
        res = c(notCygdrive, res)[1]
      }
      res
    }
    subfun = function(xx) {
      gsub("^([[:alpha:]])(:)", '/cygdrive/\\1', xx)
    }
  } else {
    whichfun = function(xx) 
    gsub("xetex", "xelatex", normalizePath(Sys.which(xx),'/'))    
    subfun = function(xx) xx
  }
  
  if(fullpath) {
    before = list(
      REXE = subfun(file.path(R.home("bin"), "R")),
      PANDOC = normalizePath(Sys.which("pandoc"), winslash='/'),
      XELATEX = whichfun("xelatex"),
      BIBER = whichfun("biber"),
      RM = whichfun("rm"),
      CP = whichfun("cp"),
      SED = whichfun("sed"),
      PANDOCCROSSREF = whichfun('pandoc-crossref')
      )
  } else {
    before = list(REXE = 'R', PANDOC='pandoc', XELATEX = 'xelatex',
      BIBER='biber', RM='rm', CP='cp', PANDOCCROSSREF='pandoc-crossref', 
      SED='sed')
  }

  before = c(before, list(
    pandocTo = c('beamer', 'latex')[2L-beamer],
    NEWCITEPROC = '--citeproc',
    OLDCITEPROC = '--filter=pandoc-citeproc',
    CITEPROC = '$(NEWCITEPROC)'
    ))

   # if can find pandoc-crossref, create the filter string
  if(nchar(Sys.which(before$PANDOCCROSSREF))) {
    before$PANDOCCROSSREF = paste("--filter=", before$PANDOCCROSSREF, sep="")
  } else {
      # if can't find it, don't do pandoc-crossref
    before$PANDOCCROSSREF = ''
  }

  before = paste(paste(names(before), before, sep='='), collapse='\n')
  before = paste(before, '\n\nall: ', paste(x, collapse=' '), '\n\n')

  makefileTemplate = scan(
    system.file('extsrc','knitrMakefile', package='Pmisc'),
    what='a', sep='\n', quiet=TRUE)

  makefileTemplate = makefileTemplate[
  seq(grep("## end of header", makefileTemplate)+1, 
    length(makefileTemplate))
  ]

  makefileTemplate = paste("\n", makefileTemplate, sep='')
  makefileTemplate = gsub("^\n\t", "\t", makefileTemplate)
  makefileTemplate = gsub("^\n[#]", "## ", makefileTemplate)


  result = paste(
    c(before, makefileTemplate), collapse='\n')


#base::library('base');.libPaths(file.path(Sys.getenv('USERPROFILE'),'Documents','R','win-library','3.5'));
#"base::library('base');.libPaths(.Library);"


  if(is.function(output)) {
    output(result)
  }
  if(is.character(output)) {
    cat(result, file=output)
  }
  invisible(result)
}