
legendBreaks = function(pos, breaks, ...){
	ldots = list(...)
	defaults = list(pch=15,  x=pos,bg="white",
		inset=0.001,cex=1)
	
	if(is.list(breaks)){
		if(length(breaks$legendCol) ) {
			ldots$col=breaks$legendCol
		}
		if( length(breaks$breaks)) {
			breaks=breaks$breaks
		}
	}

	if(length(breaks)-1 == length(ldots$col)) {
		defaults$adj=c(0,0.5)
		ldots$col = c(NA, ldots$col)
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
	
	
	
	do.call(legend, ldots)
	
}
