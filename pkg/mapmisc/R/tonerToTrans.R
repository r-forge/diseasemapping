
tonerToTrans = function(x, power=0.5, col='black') {
	if(!any(names(x)=='stamen.tonerRed'))
		warning("x doesn't seem to be stamen-toner")

	newTrans = ((255-values(x[['stamen.tonerRed']]))/255)^power
	newTrans[is.na(newTrans)] = 0
	
	col = col2rgb(col)[,1]
	
	newCol = grDevices::rgb(
			col['red'],
			col['green'],
			col['blue'],
			newTrans
	)
	newCol = factor(newCol)
	levels(newCol) = gsub("FFFFFFFF", "FFFFFF00", levels(newCol))
	result = raster(x)
	names(result) = 'stamen.toner'
	values(result) = as.numeric(newCol)
	result@legend@colortable = c(NA,levels(newCol))
	result
}
