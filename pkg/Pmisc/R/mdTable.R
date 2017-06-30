#' Markdown tables
#' 
#' @description Markdown tables with pandoc-crossref support
#'
#' @param x a matrix or data frame
#' @param col.just horizontal justification of text in cells
#' @param guessGroup attempt to guess cgroup and rgroup, for reshape2 output
#' @param ... other arguments for \code{\link[Hmisc]{latex}}, \code{\link[htmlTable]{htmlTable}} or \code{\link[knitr]{kable}}.
#' @param mdToTex if \code{FALSE} convert to markdown table, return \code{x} otherwise
#' @details \code{\link[Hmisc]{latex}} is called if \code{mdToTex} is \code{TRUE}.  Otherwise, 
#' \code{\link[htmlTable]{htmlTable}} is called if arguments such as \code{cgroup} are present
#' and \code{\link[knitr]{kable}} is called to produce a markdown table.
#' 
#' Captions with labels suitable for pandoc-crossref are added.  If not provided, the
#' label is \code{tbl:} followed by the chunk label.
#' 
#' When \code{mdToTex} is \code{'auto'} (the default), it is set to
#' \code{any(commandArgs()=='mdToTex', na.rm=TRUE)}
#' @return character string which knitr prints 'asis'
#' @examples
#' mytable = as.data.frame(matrix(runif(6),,3))
#' names(mytable) = month.name[1:ncol(mytable)]
#' 
#' cat(mdTable(x=mytable, digits=3,
#'   caption = 'the table', caption.loc='bottom',
#' mdToTex=FALSE))
#' 
#' @export
mdTable = function(x, col.just = 'r', guessGroup=FALSE, ..., mdToTex = 'auto') {
  
  if(identical(mdToTex, 'auto'))
    mdToTex = any(commandArgs()=='mdToTex', na.rm=TRUE)
  
  dots = c(
      list(x=x, 
          col.just = rep_len(col.just, ncol(x))),
      list(...))
  
  if(guessGroup) {
    
    rTable= table(x[,1])
    xFirstRowUnique = as.character(x[!duplicated(x[,1]), 1])
    if(length(xFirstRowUnique) != length(names(rTable))) {
      warning("first column doesn't appear to be a grouping")
    } 
    if( (!all(rTable==1)) & all(
        xFirstRowUnique %in% names(rTable)
    ) ) {
      rTable = drop(as.matrix(rTable))[xFirstRowUnique]
      if(!length(dots$rgroup)) {
        dots$rgroup = names(rTable)
      }
      if(!length(dots$n.rgroup)) {
        dots$n.rgroup = rTable
      }
      if(!length(dots$rowname)) {
        dots$rowname = x[,2]
      }
      rownames(dots$x) = NULL
      dots$x = dots$x[,seq(3, ncol(dots$x))]
    }
    
    # if first column is unique values and character or factor, override row names
    if(all(rTable==1) & !length(dots$rowname)) {
      if(is.character(dots$x[,1]) | is.factor(dots$x[,1])) {
        dots$rowname = as.character(dots$x[,1])
        dots$x= dots$x[,-1]
      }
    }
    
    # look for row and column groups
    cTable= table(gsub("_([[:alnum:]]|[[:space:]]|[.])+$", "", colnames(dots$x)))
    if( !all(cTable==1) ) {
      if(!length(dots$cgroup)) {
        dots$cgroup = names(cTable)
      }
      if(!length(dots$n.cgroup)) {
        dots$n.cgroup = cTable
      }
      if(!length(dots$colheads)) {
        dots$colheads = gsub(
            paste(
                paste(names(cTable), '_?', sep=''), 
                collapse='|'), "", 
            colnames(dots$x))
      }
    }
    
  } # end guessing groups
  
  
  theLabel = c(dots$label, 
      paste("tbl:", 
          c(knitr::opts_current$get()$label, 'labelMissing')[1],
          sep='')
  )[1]
  
  missingNames = which(!nchar(names(dots)))
  
  if(identical(mdToTex, TRUE)) {
    # produce latex table using Hmisc::latex
    requireNamespace("Hmisc", quietly=TRUE)
    newNames = names(formals(Hmisc::latex))[missingNames]
    names(dots)[missingNames] = newNames
    
    if( (!'object' %in% names(dots)) & ('x' %in% names(dots)) ) {
      names(dots) = gsub("^x$", "object", names(dots))
    }
    
    dots$label = theLabel
    dots$file = ''
    if(!length(dots$title)) dots$title = ''
    
    
    res = utils::capture.output(invisible(
            do.call(Hmisc::latex, dots)))
    res = res[grep("^%", res, invert=TRUE)]
    if(!identical(dots$na.blank, FALSE))
      res = gsub("(&)[[:space:]]*NA[[:space:]]*(&|\\\\)", "\\1  \\2", res)
    
  } else { # not mdToTex
    # produce HTML table with htmlTable::htmlTable or knitr::kable
    
    # some options will be ignored if knitr::kable is used
    getRidForKable = c(
        'caption.loc', 'caption.lot', 'pos.caption', 
        'label', 'row.label', 'title', 'fig.pos',
        'table.env', 'center',
        'booktabs','ctable', 'where')
    
    
    # use kable if there aren't rgroup and cgroup commands
    if(all(
        names(dots) %in% 
            c(names(formals(knitr::kable)), getRidForKable) 
    ) ) {
      dots$format = 'markdown'
      res = as.character(do.call(
              knitr::kable, 
              dots[setdiff(names(dots), getRidForKable)]
          ))
      res = c(res, '', paste(
              ": ", dots$caption, 
              " {#", theLabel, "}\n\n", sep='')) 
    } else {
      # use htmlTable
#      newNames = names(formals(htmlTable::htmlTable))[missingNames]
      newNames = names(formals(utils::getFromNamespace('htmlTable.default', 'htmlTable')))
      names(dots)[missingNames] = newNames[missingNames] 
      
      formatArgs = intersect(
          names(formals(format.default)),
          names(dots))
      if(length(formatArgs)>1) { # more than x
        dots$x =  do.call(format, dots[formatArgs])
        dots= dots[c('x', setdiff(names(dots), formatArgs))]
      }
      # remove leading white space, which causes problems for html tables
      dots$x = apply(dots$x, 2, trimws)
      
      
      dots$label = theLabel
      dots$file = ''
      
      
      
      # remove the caption, it will be added later
      theCaption = dots$caption
      dots = dots[names(dots) != 'caption']
      
      
      # names to convert from hmisc names to htmlTable names
      convertToHtmlTable = c(align = 'col.just', 
          header = 'colheads', 
          rnames = 'rowname')
      for(D in names(convertToHtmlTable)) {
        dots[[D]] = dots[[ convertToHtmlTable[D] ]]
      }
      
      res = do.call(htmlTable::htmlTable, dots)
      res = unlist(strsplit(res, '\n'))
#
#        capPos = min(grep("<table", res))+1
#      } else {
#        capPos = max(grep("[<][/]table", res))
#      }
      
      # a hack, create empty table with a caption
      
  if(identical(dots$pos.caption, 'top')) {
  res = c(
          '<div>',
          '|   |   |\n|---|---|',
          '<span style="display:inline-block; width: 50em"> <span>| \n',   
           paste(
               ': ', theCaption, ' {#tbl:',
              c(knitr::opts_current$get()$label, 'labelMissing')[1],
              '}\n', sep=''),
          res,
          '</div>\n'
      )
    } else {
      res = c(
          '<div>',
          res,
          '\n|   |   |\n|---|---|',
          '<span style="display:inline-block; width: 50em"> <span>| \n',   
          paste(
              ': ', theCaption, ' {#tbl:',
              c(knitr::opts_current$get()$label, 'labelMissing')[1],
              '}\n', sep=''),
          '</div>\n'
      )
      
    }
      
#      res[capPos] = paste("<caption>", theCaption, "</caption>\n",
#          res[capPos])
    } # end use htmltable  
  } # end not latex
  res = paste(res, '\n', sep='')
  knitr::asis_output(res)
}
