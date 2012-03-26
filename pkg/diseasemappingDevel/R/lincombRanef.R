lincombRanef = function(x, ranefName = "cellID") {
	
	if(is.vector(x)) 
		x = as.data.frame(x)
	
	if(! ranefName %in% names(x))
	 x[,ranefName] = seq(1, nrow(x))

 
 result = list()
	
	
	for(D in 1:nrow(x)){
		result[[D]]=list()
		result[[D]][[1]] = list(
				'(Intercept)' = list(weight=1)
				
				)
		for(Dcol in 1:ncol(x)) {
			result[[D]][[Dcol+1]] = list(
					list(weight=
							x[D,Dcol])
					)
			names(result[[D]][[Dcol+1]]) = colnames(x)[Dcol]		
		}
		N = length(result[[D]])
		
		result[[D]][[N]][[ranefName]]$idx = result[[D]][[N]][[ranefName]]$weight 				
		result[[D]][[N]][[ranefName]]$weight = 1
	}
	result
	
}