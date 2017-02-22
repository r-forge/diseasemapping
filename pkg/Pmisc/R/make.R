#' @export
make = function(x, 
    suffix=NULL, beamer=FALSE,
    run=FALSE, ...) {
  
  if(!is.null(suffix))
    x = paste(
        tools::file_path_sans_ext(basename(x)),
        suffix, sep='.')
  
  theString= paste(' -f', 
      system.file('src','knitrMakefile', package='Pmisc'),
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
Makefile = function(x, suffix = NULL,
    beamer = FALSE, output='Makefile', cygwin=FALSE) {
  
  if(!is.null(suffix))
    x = paste(
        tools::file_path_sans_ext(basename(x)),
        suffix, sep='.')
  
  if(cygwin){
    whichfun = function(x) {
      res = system(paste(
              Sys.which("bash"),
              " -c 'which --all ", x, "'",
              sep=''
          ), intern=TRUE)
      if(length(res > 1)) {
          notCygdrive = grep("cygdrive", res, value=TRUE, invert=TRUE)
          res = c(notCygdrive, res)[1]
      }
      res
    }
    subfun = function(x) {
      gsub("^([[:alpha:]])(:)", '/cygdrive/\\1', x)
    }
  } else {
    whichfun = function(x) Sys.which(x)
    subfun = function(x) x
  }
  
  before = list(
      REXE = subfun(file.path(R.home("bin"), "R")),
      PANDOC = normalizePath(Sys.which("pandoc"), winslash='/'),
      XELATEX = whichfun("xelatex"),
      BIBER = whichfun("biber"),
      RM = whichfun("rm"),
      pandocTo = c('beamer', 'latex')[2L-beamer],
      DOCXTEMPLATE = 
          system.file('src/template.docx', package='Pmisc'),
      ODTTEMPLATE = 
          system.file('src/template.odt', package='Pmisc')
  )
  
  before = paste(paste(names(before), before, sep='='), collapse='\n')
  before = paste(before, '\n\nall: ', paste(x, collapse=' '), '\n\n')
  
  makefileTemplate = scan(
      system.file('src','knitrMakefile', package='Pmisc'),
      what='a', sep='\n', quiet=TRUE)
  
  makefileTemplate = makefileTemplate[
      seq(min(grep("^[.]PRECIOUS[:]", makefileTemplate)), 
          length(makefileTemplate))
  ]

  makefileTemplate = paste("\n", makefileTemplate, sep='')
  makefileTemplate = gsub("^\n\t", "\t", makefileTemplate)
  makefileTemplate = gsub("^\n[#]", "## ", makefileTemplate)
  
  # template variables
  makefileTemplate = gsub("[$][$][(][$][(]((ODT|DOCX)TEMPLATE)[)][)]", "$(\\1)",
       makefileTemplate)
  
  result = paste(
      c(before, makefileTemplate), collapse='\n')
  
  
  '.PRECIOUS: %.md %.tex
      
      %.md: %.Rmd
      $(REXE) -e "knitr::knit(\'$<\', encoding=\'UTF-8\')" $(Rargs)
      
      %.tex: %.md
      $(PANDOC) --standalone --smart --biblatex $(pandocArgs) --to=$(pandocTo) --output=$@ $<
      
      %.bcf: %.tex
      $(XELATEX) $<
      
      %.bbl: %.bcf
      $(BIBER) $<	
      
      %.pdf: %.tex %.bbl
      $(XELATEX) -interaction=nonstopmode $<;
      $(XELATEX) -interaction=nonstopmode $<
      
      %.html: %.md
      $(PANDOC) --smart --standalone --mathjax --filter=pandoc-citeproc $(pandocArgs) --to=html5 --output=$@ $<
      
      %.rtf: %.md
      $(PANDOC) --smart --standalone --filter=pandoc-citeproc $(pandocArgs) --output=$@ $<
      
      %.odt: %.md
      $(PANDOC) --smart --standalone --filter=pandoc-citeproc $(pandocArgs) --reference-docx=$(odtTemplate) --output=$@ $<
      
      %.docx: %.md
      $(PANDOC) --smart --standalone --filter=pandoc-citeproc $(pandocArgs) --reference-docx=$(docxTemplate) --output=$@ $<
      
      clean:
      $(RM) *.run.xml* *.blg *.out *.log *.aux *.bcf *.bbl *.nav *.toc *.vrb'
  
  
  if(is.function(output)) {
    output(result)
  }
  if(is.character(output)) {
    cat(result, file=output)
  }
  invisible(result)
}