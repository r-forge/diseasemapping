map.new = function(x) {
	

	
	xpoints = t(bbox(extent(x)))

	xpoints = SpatialPoints(xpoints)
	thecrs = try(proj4string(x), silent=TRUE) 
	if(class(thecrs)!="try-error")
		proj4string(xpoints) = CRS(thecrs)

	
	par(mar=c(0,0,0,0))
	plot(xpoints,pch=NA)

	
	return(invisible())
		
}

legendBreaks = function(pos, col, breaks, ...){
	ldots = list(...)
	defaults = list(pch=15, adj=c(0,0.5), x=pos,
			inset=0.001,cex=1)
	for(D in names(defaults))
		if(is.null(ldots[[D]]))
			ldots[[D]] = defaults[[D]]

	defaults = list(pt.cex=2.5*ldots[["cex"]])
	for(D in names(defaults))
		if(is.null(ldots[[D]]))
			ldots[[D]] = defaults[[D]]
	
	
	ldots$legend = rev(breaks)
	ldots$col = c(rev(col),NA)

	
	# get rid of transparency in col
	withTrans = grep("^#[[:xdigit:]]{8}$", ldots$col)
	ldots$col[withTrans] = gsub("[[:xdigit:]]{2}$", "", ldots$col)
	
 
	
	do.call(legend, ldots)
}
