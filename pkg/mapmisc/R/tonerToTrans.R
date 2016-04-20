
tonerToTrans = function(x, power=0.5) {
	if(!any(names(x)=='stamen.tonerRed'))
		warning("x doesn't seem to be stamen-toner")
	
	newCol = grDevices::rgb(
			0,
			0,
			0,
			((255-values(x[['stamen.tonerRed']]))/255)^power
	)
	newCol = factor(newCol)
	levels(newCol) = gsub("FFFFFFFF", "FFFFFF00", levels(newCol))
	result = raster(x)
	names(result) = 'stamen.toner'
	values(result) = as.numeric(newCol)
	result@legend@colortable = c(NA,levels(newCol))
	result
}
