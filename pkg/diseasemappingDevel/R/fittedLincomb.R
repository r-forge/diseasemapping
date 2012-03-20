fittedLincomb = function(inlaResult){
	
	allterms = attributes(terms(inlaResult$formula))$term.labels
	theterms = allterms[-grep("^f\\(", allterms)]
	theRE = names(inlaResult$summary.random)
	if(length(theRE)>1)
		warning("more than one random effect", theRE)
	
	theUnique = which(!duplicated(inlaResult[[theRE]]))
	forCellID = matrix(NA, dim(inlaResult$summary.random[[theRE]])[1],length(theUnique))
	
	for(D in seq(1, length(theUnique))) {
    	forCellID[inlaResult$data[theUnique[D],theRE],D] = 1
	}
	
	
	forLincomb = list(cellID=forCellID, "(Intercept)"=rep(1,length(theUnique)))
	names(forLincomb)[1] = theRE
	
	theLincomb = inla.make.lincombs(inla.uncbind(as.matrix(
							inlaResult$data[theUnique,theterms])), 
			forLincomb)
	
	theLincomb
	
}
