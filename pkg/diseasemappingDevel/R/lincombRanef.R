lincombRanef = function(x, ranefName = "cellID", intercept=T) {
	
	if(length(x)==1) {
		x = 1:x
	}
	
	if(is.vector(x)) {
		x = as.data.frame(x)
		names(x) = ranefName
	}
		
	if(! ranefName %in% names(x))
	 x[,ranefName] = seq(1, nrow(x))

 
 result = list()
	
	
	for(D in 1:nrow(x)){
		result[[D]]=list()
		for(Dcol in 1:ncol(x)) {
			result[[D]][[Dcol]] = list(
					list(weight=
							x[D,Dcol])
					)
			names(result[[D]][[Dcol]]) = colnames(x)[Dcol]		
		}
	
		

		N = length(result[[D]])
		
		result[[D]][[N]][[ranefName]]$idx = result[[D]][[N]][[ranefName]]$weight 				
		result[[D]][[N]][[ranefName]]$weight = 1
		
		
		if(intercept) {
			result[[D]] = c(result[[D]],
				list(list(
					'(Intercept)' = list(weight=1)
			
					) 
				))
		}
		
	}
	result
	
}