
legendBreaks = function(pos,
    breaks,
    col=breaks$col,
    legend=breaks$breaks,
    rev=TRUE,
    outer=TRUE,
    pch=15,
    cex=par('cex'),
    pt.cex=2.5*cex,
    inset=0.05,
    text.col=par('fg'),
    title.col=text.col,
    adj=0,
    ...){
  
  if(rev){
    col=rev(col)
    legend=rev(legend)
  }
  
  if(length(col) == (length(legend)-1)) {
    col = c(NA, col)
  }
  
  pch = c(NA,
      pch[round(seq(1, length(pch), len=length(legend)-1))]
  )
  
  # get rid of transparency in col
  withTrans = grep("^#[[:xdigit:]]{8}$", col)
  col[withTrans] = gsub("[[:xdigit:]]{2}$", "", col)

  if(outer){
    oldxpd = par("xpd")
    par(xpd=NA)
    fromEdge = matrix(par("plt"), 2, 2, 
        dimnames=list(c("min","max"), c("x","y")))
    propIn = apply(fromEdge, 2, diff)
    if(is.character(pos)) {
      forInset = c(0,0)
      if(length(grep("left$", pos))){
        forInset[1] = -fromEdge["min","x"]					
      } else if(length(grep("right$", pos))){
        forInset[1] = fromEdge["max","x"]-1					
      }
      if(length(grep("^top", pos))){
        forInset[2] = -fromEdge["min","y"]					
      } else if(length(grep("^bottom", pos))){
        forInset[2] = fromEdge["max","y"]-1					
      }
      
      inset = forInset/propIn + inset
    }
  }
  
  result=legend(
      pos,
      legend=legend,
      col=col,
      pch=pch,
      pt.cex=pt.cex,
      inset=inset,
      cex=cex,
      text.col='#FFFFFF00',
      ...
      )
      
      diffy = diff(result$text$y)/2
      diffy = c(diffy,diffy[length(diffy)])
      result$text$y = result$text$y + diffy
      
      if(par("xlog")) result$text$x = 10^result$text$x
      if(par("ylog")) result$text$y = 10^result$text$y
      
      
      text(result$text$x, result$text$y,
          legend, col=text.col,adj=adj)   
      
      par(xpd=oldxpd)
      
      return(invisible(result))
}

legendBreaksOld = function(pos, breaks, outer=TRUE,...){
	ldots = list(...)

	defaults = list(pch=15,  x=pos,bg="white",
		inset=0.001,cex=1)
	
if(outer){
	oldxpd = par("xpd")
	par(xpd=NA)
	fromEdge = matrix(par("plt"), 2, 2, 
			dimnames=list(c("min","max"), c("x","y")))
	propIn = apply(fromEdge, 2, diff)
	if(is.character(pos)) {
		inset = c(0,0)
#		if(length(grep("^bottom", pos))){
#			inset[2] = -fromEdge["min","y"]					
#		} else if(length(grep("^top", pos))){
#			inset[2] = fromEdge["max","y"]-1					
#		}
		
		if(length(grep("left$", pos))){
			inset[1] = -fromEdge["min","x"]					
		} else if(length(grep("right$", pos))){
			inset[1] = fromEdge["max","x"]-1					
		}
		inset = inset/propIn
		defaults$inset = inset+ 0.01
	}
} else {
	defaults$inset = 0.01
}


	if(is.list(breaks)){
		if(length(breaks$col) ) {
			ldots$col=breaks$col
		}
		if( length(breaks$breaks)) {
			breaks=breaks$breaks
		}
	}

	
	
	for(D in names(defaults))
		if(is.null(ldots[[D]]))
			ldots[[D]] = defaults[[D]]
	
	defaults = list(pt.cex=2.5*ldots[["cex"]])
	for(D in names(defaults))
		if(is.null(ldots[[D]]))
			ldots[[D]] = defaults[[D]]
	
	
	ldots$legend = rev(breaks)
	ldots$col = rev(ldots$col)
	
	
	# get rid of transparency in col
	withTrans = grep("^#[[:xdigit:]]{8}$", ldots$col)
	ldots$col[withTrans] = gsub("[[:xdigit:]]{2}$", "", ldots$col)
	
	if(length(breaks)-1 == length(ldots$col)) {
		ltext = format(ldots$legend,
				width=max(nchar(ldots$legend)))
		ldots$col = c(NA, ldots$col)
		ldots$legend = rep(NA, length(ldots$col))
		ldots$text.width = max(strwidth(ltext))
	} 
	
	result=do.call(legend, ldots)
 
	
	if(all(is.na(ldots$legend))) {
		x=ldots
		if(!is.null(x$title))
			x$title = NA
		x$legend=ltext
		x$x = result$rect$left
		x$y = result$rect$top - 0.75*par("cxy")[2]
		x$adj = c(-0.2, 0) 
		x$col=NA
		x$bty='n'
 
		do.call(legend,x)
	}

	par(xpd=oldxpd)
	
	return(invisible(result))
}
