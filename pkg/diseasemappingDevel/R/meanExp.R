meanExp = function(marginals) {
	
	result = meanNatural = rep(NA, length(marginals))
	
	for(D in 1:length(marginals)) {
		
		theX = marginals[[D]][,1]
		
		binwidth = apply(matrix(c(0,rep(diff(theX)/2,rep(2,length(theX)-1)),0),nrow=2), 2, sum)
		
		result[D] = sum(binwidth*exp(theX) * marginals[[D]][,2])
		
		
	}
	names(result) 	=  names(marginals)
	return(result)
}