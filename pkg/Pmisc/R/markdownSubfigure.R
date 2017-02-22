#' @export
hook_plot_mdsubfig = function(x, options) {
# make multiple plots as pseudo-subfigures, in a table
  
  
  if (options$fig.show == 'animate') return(knitr::hook_plot_html(x, options))
  
  base = knitr::opts_knit$get('base.url') 
  if(is.null(base)) base=''
  cap = options$fig.cap
  scap = options$fig.subcap
  if(is.null(cap)) cap=''
  if(is.null(scap)){
    scap = cap
  }
  
  result = 
      sprintf('![%s](%s%s) ', scap, base, x)#knitr:::.upload.url(x))    
  
  if(any(options$fig.ncol==0)){
    return(result)
  }
  
  
  fig.num =options$fig.num
  if(is.null(fig.num))
    fig.num=1
  
  
  
  
  if(length(scap)) {
    fig.cur =options$fig.cur
    if(is.null(fig.cur))
      fig.cur=1
    fig.ncol = options$fig.ncol
    if(is.null(fig.ncol)) {
      fig.ncol = 1
    }
    
    Drow = floor((fig.cur-1)/fig.ncol)+1
    Dcol = fig.cur - (Drow-1) * fig.ncol 
    
    
    if( (Dcol == fig.ncol ) | (fig.cur == fig.num) ) { 
      # we're at the end of a column or the last plot
      Dend = fig.ncol*Drow

 
      
      result = paste(result, 
          "\n\n", sep=""
      )
      
    } else {
      result = paste( result, "\n", sep="")  	
    }
    
    if(fig.cur==1) {
      result = paste('<div id="',
          options$fig.lp, 
          options$label, 
          '">\n',
          result, sep="")
    }
    
    if(fig.cur == fig.num)
      result = paste(result, cap, "\n</div>\n\n",
          sep="")	
    
  }
  result
  
}
