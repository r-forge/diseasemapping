fittedLincomb = function(inlaResult){
	
	allterms = attributes(terms(inlaResult$formula))$term.labels
	theterms = allterms[-grep("^f\\(", allterms)]
	theRE = names(inlaResult$summary.random)
	if(length(theRE)!=1)
		warning("more than one random effect", theRE)
	
	theUnique = which(!duplicated(inlaResult$data[[theRE]]))
	forCellID = matrix(NA, dim(inlaResult$summary.random[[theRE]])[1],length(theUnique))
	
	
	for(D in seq(1, length(theUnique))) {
    	forCellID[inlaResult$data[theUnique[D],theRE],D] = 1
	}
	
	
	forLincomb = list(cellID=forCellID, "(Intercept)"=rep(1,length(theUnique)))
	names(forLincomb)[1] = theRE
	
	#return(list(theUnique, theterms, forLincomb,inla.uncbind(as.matrix(
	#								inlaResult$data[theUnique,theterms]))))
	
	theMatrix = as.matrix(inlaResult$data[theUnique,theterms])
	theMatrix2 = inla.uncbind(theMatrix)
	
	assign("theMatrix2QQ", theMatrix2,pos=1)
	
	
	theLincomb = inla.make.lincombs(theMatrix2QQ, forLincomb)
	
	rm(theMatrix2QQ, pos=1)
	
	theLincomb
	
}
