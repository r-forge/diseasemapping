
legendBreaks = function(pos, col, breaks, ...){
	ldots = list(...)
	defaults = list(pch=15,  x=pos,
			inset=0.001,cex=1)
	
	if(length(breaks)-1 == length(col)) {
		defaults$adj=c(0,0.5)
		col = c(NA, col)
	} 
	
	
	for(D in names(defaults))
		if(is.null(ldots[[D]]))
			ldots[[D]] = defaults[[D]]
	
	defaults = list(pt.cex=2.5*ldots[["cex"]])
	for(D in names(defaults))
		if(is.null(ldots[[D]]))
			ldots[[D]] = defaults[[D]]
	
	
	ldots$legend = rev(breaks)
	ldots$col = rev(col)
	
	
	# get rid of transparency in col
	withTrans = grep("^#[[:xdigit:]]{8}$", ldots$col)
	ldots$col[withTrans] = gsub("[[:xdigit:]]{2}$", "", ldots$col)
	
	
	
	do.call(legend, ldots)
	
}
