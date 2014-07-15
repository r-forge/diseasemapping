
legendBreaks = function(pos, breaks, outer=TRUE,...){
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
