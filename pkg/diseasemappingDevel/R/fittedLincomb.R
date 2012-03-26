fittedLincomb = function(inlaResult, fixedValues = NULL){
	
	allterms = attributes(terms(inlaResult$formula))$term.labels
	theterms = allterms[-grep("^f\\(", allterms)]
	
	
	if(!is.null(fixedValues))
		theterms = theterms[!theterms %in% names(fixedValues)]
		
	theRE = names(inlaResult$summary.random)
	if(length(theRE)!=1)
		warning("more than one random effect", theRE)
	
	theUnique = which(!duplicated(inlaResult$data[[theRE]]))

	

	
	for(D in names(fixedValues)) {
		if(length(fixedValues[[D]])==1) {
			forLincomb[[D]] = rep(fixedValues[[D]], length(theUnique))	
		} else {
			forLincomb[[D]] = fixedValues[[D]]
		}
	}
	
	
	
	
	#return(list(theUnique, theterms, forLincomb,inla.uncbind(as.matrix(
	#								inlaResult$data[theUnique,theterms]))))
	
	theDF = inlaResult$data[theUnique,theterms]
	theDF = theDF[,lapply(theDF, class)=="numeric"]
	theMatrix = as.matrix(theDF)
	theMatrix = cbind(theMatrix, "(Intercept)" = 1)
	theMatrix2 = inla.uncbind(theMatrix)
	
	assign("theMatrix2QQ", theMatrix2,pos=1)
	

	theLincomb = inla.make.lincombs(theMatrix2QQ)

	idVec = inlaResult$data[theUnique,theRE]
	
	for(Dcell in 1:length(theLincomb)) {
	theLincomb[[Dcell]] = c(theLincomb[[Dcell]], 
    		list(list(list(idx=idVec[Dcell], weight=1)))
			)
	Dre = length(theLincomb[[Dcell]])
			
	names(theLincomb[[Dcell]][[Dre]])=theRE		
	
	}
 
	attr(theLincomb,"idx") = inlaResult$data[[theRE]][theUnique]
 
	return(theLincomb)
	
}
