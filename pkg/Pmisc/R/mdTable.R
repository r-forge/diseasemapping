#' Markdown tables
#' 
#' @description Markdown tables with pandoc-crossref support
#'
#' @param x a matrix or data frame
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
mdTable = function(x, ..., mdToTex = 'auto') {
  
  if(identical(mdToTex, 'auto'))
    mdToTex = any(commandArgs()=='mdToTex', na.rm=TRUE)
  
  dots = c(list(x=x), list(...))
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
      newNames = names(formals(htmlTable::htmlTable))[missingNames]
      names(dots)[missingNames] = newNames
      
      formatArgs = intersect(
        names(formals(format.default)),
        names(dots))
      if(length(formatArgs)>1) { # more than x
        dots$x = do.call(format.default, dots[formatArgs])
        dots= dots[c('x', setdiff(names(dots), formatArgs))]
      }
      
      dots$label = theLabel
      dots$file = ''
      
      # remove the caption, it will be added later
      theCaption = dots$caption
      dots = dots[names(dots) != 'caption']
      
      res = do.call(htmlTable::htmlTable, dots)
      res = unlist(strsplit(res, '\n'))
      if(identical(dots$pos.caption, 'top')) {
        capPos = min(grep("<table", res))+1
      } else {
        capPos = max(grep("[<][/]table", res))
      }
      res[capPos] = paste("<caption>", theCaption, "</caption>\n",
        res[capPos])
    } # end use htmltable  
  } # end not latex
  res = paste(res, '\n', sep='')
  knitr::asis_output(res)
}
