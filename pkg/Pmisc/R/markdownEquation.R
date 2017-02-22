#' Markdown equation array
#' 
#' @description Multiline equation in a markdown table
#'
#' @param x a string with latex code for the equation
#' @param mdToTex if \code{FALSE} convert to markdown table, return \code{x} otherwise
#' @details A multi-line equation in an eqnarray or align document is converted to markdown code suitable for conversion to docx.
#' @examples
#' cat(eqnarray("
#' \begin{align*}
#' Y_i \sim & N(\mu_i, \tau^2)\\
#' \mu_i= & X_i \beta
#' \end{align*}
#' ", FALSE))
#' 
#' @export
eqnarray = function(x, mdToTex=FALSE) {
  
  if(mdToTex) {
    res = x 
  } else {
    
    res = gsub(
        "[$][$]|([[:punct:]](begin|end)[[:punct:]](eqnarray|align(ed)?)[*]?[[:punct:]])", 
        "", x)
    res = gsub("[[:space:]]*&[[:space:]]*", " ", res)
    res = gsub(" *[\\][\\] *\n? *", "$ |  |\n| $", res)
    
    res = gsub("\n *[|] *[$][[:space:]]*$", "\n", res)
    
    res = gsub("^[[:space:]]*([|] *[$]){0}", "| $", res)
    res = gsub("([$] *[|] *[|])?[[:space:]]*$", "$ |  |\n\n", res)
    
    res = paste(
        '|      |       |\n|:-----:|:------|\n', 
        gsub("^([[:space:]]+)?|[|]$", "", res), 
        sep='')
  }
  knitr::asis_output(res)
}
