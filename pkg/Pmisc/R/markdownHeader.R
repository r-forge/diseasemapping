#' @export
subcaptionCommands = c(
  '\\expandafter\\def\\csname ver@subfig.sty\\endcsname{}',
  '\\usepackage{subcaption}',
  '\\makeatletter\\@ifundefined{subfloat}{\\newcommand{\\subfloat}[2][need a sub-caption]{\\subcaptionbox{#1}{#2}\\, }}{ \\renewcommand{\\subfloat}[2][need a sub-caption]{ \\subcaptionbox{#1}{#2}\\, }}\\makeatother'
)


#' @export
mathCommands =c(
  '\\DeclareMathOperator{\\var}{var}',
  '\\DeclareMathOperator{\\E}{E}',
  '\\DeclareMathOperator{\\N}{N}',
  '\\DeclareMathOperator{\\MVN}{MVN}',
  '\\newcommand{\\comment}[1]{}',
  '\\newcommand{\\simiid}{\\overset{\\text{ind}}{\\sim}}',
  '\\newcommand{\\yv}{\\mathbf{y}}', 
  '\\newcommand{\\Yv}{\\mathbf{Y}}', 
  '\\newcommand{\\dens}{\\symbf{\\pi}}',
  '\\newcommand{\\tp}{^\\text{T}}'
) 


mathmlCommands = '<script type="text/x-mathjax-config">MathJax.Hub.Config({config: ["MMLorHTML.js"],jax: ["input/TeX","input/MathML","output/HTML-CSS","output/NativeMML"],extensions: ["tex2jax.js","mml2jax.js","MathMenu.js","MathZoom.js"],TeX: {extensions: ["AMSmath.js","AMSsymbols.js","noErrors.js","noUndefined.js"]}});</script><script type="text/javascript" src="https://cdn.mathjax.org/mathjax/2.0-latest/MathJax.js?config=TeX-AMS-MML_HTMLorMML"></script>'

#' @export
today = function(...){
  dots = list(...)
  if(!any(names(dots)=='format'))
    dots$format='%A %e %B %Y'
  do.call(format, c(list(x=Sys.time()), dots))
}

slideHeaderIncludes = 
#  knitr::asis_output(
#    knitr::combine_words(
        paste0(
          '\\newcommand{\\',
          c(
            'newcol}[1][0.5]{\\column{#1\\textwidth}}',
            'bcol}[1][0.5]{\\begin{columns}\\column{#1\\textwidth}}',
            'ecol}{\\end{columns}}',
            'newcolb}[2][0.5]{\\end{block}\\column{#1\\textwidth}\\begin{block}{#2}}',
            'bcolb}[2][0.5]{\\begin{columns}\\column{#1\\textwidth}\\begin{block}{#2}}',
            'ecolb}{\\end{block}\\end{columns}}')#,
#          sep=''),
#      sep='\n', and='')
  )


#' @export
markdownHeader = function(
  title=NULL, 
  author=NULL,
  date=Pmisc::today(),
  biblatex=1,
  bibliotitle = 'References',
  bibliostyle = 'authoryear,backend=biber',
  biblatexoptions = c(
    maxbibnames=20,
    maxcitenames=2,doi='false',
    isbn='false',	giveninits='true',
    uniquelist='false'),
  numbersections=TRUE,
  beamer=FALSE,
  mathCommands=FALSE,
  subcaptionCommands = !beamer,
  crossrefYaml = system.file(
    'src','crossref.yaml', package='Pmisc'),
  csl = system.file(
    'src','apa-no-doi-no-issue.csl', package='Pmisc'),
  css = system.file(
    'src','article.css', package='Pmisc'),
  ...
) {
  
  startAndEnd='---'  
#  extraHeaderIncludes = '\\@ifpackageloaded{biblatex}{\\renewbibmacro{in:}{}}{}'
  extraHeaderIncludes = NULL
  if(mathCommands)
    extraHeaderIncludes = c(Pmisc::mathCommands, extraHeaderIncludes)
  if(subcaptionCommands)
    extraHeaderIncludes = c(Pmisc::subcaptionCommands, extraHeaderIncludes)

 
  

  mode(biblatexoptions) = 'character'
  
  result = list(
    startAndEnd,
    paste('title:', title),
    paste('author:', author),
    paste('date:', date)
  )
  
  dots = list(
    biblatex=biblatex,
    'biblio-title'=bibliotitle,
    'biblio-style'=bibliostyle,
    biblatexoptions = biblatexoptions,
    numbersections=numbersections,
    numberSections=numbersections,
    graphics=TRUE,
    subfigGrid= TRUE,
    crossrefYaml = crossrefYaml,
    csl = csl, css=css,
    ...)	
  
  dots = dots[
    unlist(lapply(dots, length))>0
  ]
  
  names(dots) = gsub(
    "header[[:punct:]]?includes",
    'header-includes', names(dots),
    ignore.case=TRUE)
  
  if(length(extraHeaderIncludes)) {
    dots[['header-includes']] = c(extraHeaderIncludes,  dots[['header-includes']] )
  }

  dots[['header-includes']] = paste0(
      "  - '`",
      c(dots[['header-includes']], slideHeaderIncludes[beamer]),
      "`{=latex}'")
  dots[['header-includes']] = c(
    dots[['header-includes']],
    paste0(
      "  - '`",
      mathmlCommands,
      "`{=html}'")
    )

  
  for(D in setdiff(names(dots), 'header-includes')) {
    if(length(dots[[D]])>1) {
      dots[[D]] = knitr::combine_words(c('',
          paste(
            '  - ', 
            names(dots[[D]]),
            c("",' = ')[1+as.logical(nchar(names(dots[[D]])))],
            unlist(dots[[D]]),
            sep='')
        ), 
        sep='\n', and='')
    }
  }
  
  
  
  dotsNotHeaders = grep("^header-includes", names(dots), value=TRUE, invert=TRUE)
  
  for(D in dotsNotHeaders) {
    result[[length(result)+1]] = paste(
      D, ": ", dots[[D]], sep=''
    )	
  }
  


  # header includes
  if(beamer | any(names(dots)=='header-includes')) {

    hIncludes = dots[['header-includes']]

    toAdd <- gsub(
      '\n+', '\n', 
      knitr::combine_words(c(
          'header-includes:', 
          hIncludes),
      sep='\n', and='')
    )

    result[[length(result)+1]] = toAdd 

  }
  result[[length(result)+1]] = startAndEnd



  
  result = knitr::combine_words(c('',result,''), 
    sep='\n', and='')
  knitr::asis_output(result)  
}




