
excProb = function(marginals, threshold=0, template=NULL, templateIdCol=NULL) {
	
	if(!all(c("x","y") %in% colnames(marginals[[1]]))) {
		warning("need x and y in column names of marginals")
	}

excFunQQ = function(themat) {
	over = themat[,"x"]>threshold
	
	toInt = rbind(c(threshold, approx(themat[,"x"], themat[,"y"], threshold)$x),
			themat[over,]
			)
	
	prob = trapz(toInt[,"x"], toInt[,"y"])
	
	prob
} 

if(is.list(marginals)) {
 excProbAll = unlist(lapply(marginals, excFunQQ))
} else {
	excProbAll = excFunQQ(marginals)
}
# make sure probabilities are between zero and 1
excProbAll = pmax(0, pmin(excProbAll, 1))

if(length(grep("^Raster", class(template)))) {
	values(template) = excProbAll
	excProbAll = template
	names(excProbAll) = paste("exc", threshold, sep="")
} 

if(length(grep("(SpatialPolygonsDataFrame|SpatialPointsDataFrame)", class(template)))) {
	newcol=paste("exc", threshold, sep="")
	if(is.null(templateIdCol)) {
		template[[newcol]] = excProbAll
	} else {
	template[[newcol]] = 
			excProbAll[template@data[templateIdCol]]
	}
	excProbAll = template
} 


excProbAll
} 