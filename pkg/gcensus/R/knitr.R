#' Hooks for knitr
#' 
#' @description Hooks for including in knitr documents
#' @param x see \code{\link[knitr]{hook_pdfcrop}}
#' @param options see \code{\link[knitr]{hook_pdfcrop}}
#' 
#' @export
hook_plot_htmlsubfig = function(x, options) {
# make multiple plots as pseudo-subfigures, in a table
	
  # for debugging
  # stuff <<- list(x=x, options=options)
  
  if (options$fig.show == 'animate') return(knitr::hook_plot_html(x, options))
  
  base = knitr::opts_knit$get('base.url') 
  if(is.null(base)) base=''
  cap = options$fig.cap
  scap = options$fig.subcap
  if(is.null(cap)) cap=''
  if(is.null(scap)){
    scap = cap
  }
  
  if (is.null(w <- options$out.width) & is.null(h <- options$out.height) &
    	is.null(s <- options$out.extra) & options$fig.align == 'default') {
    result = 
    		sprintf('![%s](%s%s) ', scap, base, knitr:::.upload.url(x))
    
  } else {
    # use HTML syntax <img src=...>
    result = 
    		knitr:::.img.tag(knitr:::.upload.url(x), w, h, scap,
    				c(s, sprintf('style="%s"', knitr:::css_align(options$fig.align)))
    		)
  }
  
  if(any(options$fig.ncol==0)){
    return(result)
  }
  
  
  fig.num =options$fig.num
  if(is.null(fig.num))
    fig.num=1
  
  
  fig.subcap.all = options$params.src
  fig.subcap.all = eval(parse(text=
    					paste("list(", sub("^[[:alnum:]]+,", "", fig.subcap.all), ")")
    	))$fig.subcap
  if(length(fig.subcap.all) < fig.num)
    fig.subcap.all = c(
    		fig.subcap.all, 
    		rep(" ", fig.num - length(fig.subcap.all))
    )
  fig.subcap.all = fig.subcap.all[1:fig.num]
  
  
  if(length(fig.subcap.all)) {
    fig.cur =options$fig.cur
    if(is.null(fig.cur))
    	fig.cur=1
    fig.ncol = options$fig.ncol
    if(is.null(fig.ncol)) {
    	fig.ncol = 1
    }
    
    Drow = floor((fig.cur-1)/fig.ncol)+1
    Dcol = fig.cur - (Drow-1) * fig.ncol 
    
#	cat("\n rc ", fig.cur, Drow, " ", Dcol, "\n")
    
    if(Dcol==1) {
    	result = paste("|", result, sep="")
    }
    
    if( (Dcol == fig.ncol ) | (fig.cur == fig.num) ) { 
    	# we're at the end of a column or the last plot
    	Dend = fig.ncol*Drow
    	nextra = Dend - fig.cur
    	# pandoc doesn't like single column tables, add an extra column
    	if(fig.ncol==1) nextra = nextra + 1
    	result = paste( result,
    			paste(rep("   |  ", nextra), collapse=""),		
    			"|\n")
    	# add the subcaptions	
    	
    	Ssubcap = seq(fig.ncol*(Drow-1)+1,
    			min(c(Dend, length(fig.subcap.all)))
    	)
    	fig.subcap.all = fig.subcap.all[Ssubcap]
    	fig.subcap.all = paste(letters[Ssubcap], fig.subcap.all, sep=") ")
    	fig.subcap.all= c(	fig.subcap.all,
    			rep(" ",nextra)
    	) 
    	
    	result = paste(result, "| ",
    			paste(fig.subcap.all, 
    					collapse=" | "), 
    			" |\n", sep=""
    	)
    	
    } else {
    	result = paste( result, "|", sep="")  	
    }
    
    if(fig.cur==1) {
    	result = paste('\n\n\n|',
    			paste(rep('     ', max(c(2,fig.ncol))), collapse='|'), '|\n', 
    			"|",
    			paste(rep(':---:', max(c(2,fig.ncol))), collapse='|'),'|\n', 
    			result, sep="")
    }
    
    if(fig.cur == fig.num)
    	result = paste(result, "\n : Figure: ", cap, "\n\n\n\n", sep="")	
    
  }
  result
  
}


#' @export
hook_plot_beamer = function(x, options) {
	
  if (options$fig.show == 'animate') return(knitr::hook_plot_html(x, options))

		
  fig.ncol = options$fig.ncol
  if(is.null(fig.ncol)) {
  	fig.ncol=1
  }
	
	fig.num =options$fig.num
	if(is.null(fig.num))
		fig.num=1
	
  fig.cur =options$fig.cur
	if(is.null(fig.cur))
		fig.cur=1
	
	fig.nrow = ceiling(fig.num / fig.ncol)		
	
	Dcol = floor((fig.cur-1)/fig.nrow)+1
	Drow = fig.cur-(Dcol-1)*fig.nrow
	
#		cat(c('\n', fig.cur, fig.num, fig.nrow, fig.ncol, Drow, Dcol, '\n'))
	
	out.width  = options$out.width
	if(is.null(out.width)){
		out.width = signif(1/fig.ncol, 2)
	}
	out.width = paste("[", out.width, "]",sep='')
	
	
  base = knitr::opts_knit$get('base.url') 
	if(is.null(base)) base=''
	
  cap = options$fig.cap
  scap = options$fig.subcap
	
	scapCenter = 
  		if(is.null(scap)){
				scap = scapCenter = ''
			} else {
				scapCenter = paste("\\centering ", scap, "")
			}
	
  result = sprintf('![%s](%s%s)\\ \n\n%s\n\n', scap, base, knitr:::.upload.url(x), scapCenter)
  
	if(Drow == 1 & fig.cur > 1) {
		result = paste("\n\\newcol", out.width, "\n\n", result, sep="")
	}
	
	if(any(nchar(cap)>1)) {
		if( fig.ncol > 1) {
			if(fig.cur == 1 ) 
				result = paste("\n\\bcol", out.width, "\n\n", result, sep="")
			if(fig.cur == fig.num) 
				result = paste(result, "\n\\ecol\n\n", sep="")
		}
		if(fig.cur == 1 ) {
			result = paste( "\n", options$fig.start," ", cap,  "\n\n", result, sep="")	
		}
	}
	
  result
	
}
