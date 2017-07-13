if(FALSE) {
  # creating the precision matrix entries in maternGmrfPrec.R
  
  params=c(oneminusar=0.1, shape=0, conditionalVariance=1)
  
  Ngrid = 5
  myGrid = geostatsp::squareRaster(extent(-Ngrid,Ngrid,-Ngrid,Ngrid), 1+2*Ngrid)
  cornerSeq = list(row = seq(1,Ngrid+1), col = seq(Ngrid+1, len=Ngrid+1))
  midcell = cellFromRowCol(myGrid, 
      round(nrow(myGrid)/2), round(ncol(myGrid)/2)) 
  
  precMat =maternGmrfPrec(myGrid, param=params) 
  precMatChar = as.character(precMat)
  precMatChar[precMatChar==as.character(-1/attributes(precMat)$info$theo['a'])] = '-a'
  precMatChar = matrix(precMatChar, nrow(precMat), ncol(precMat))
  
  precForYacas = apply(precMatChar,1,paste,collapse=',')
  precForYacas = paste(precForYacas, collapse='},{')
  yacas0 = Ryacas::yacas(paste("prec0:= { {", precForYacas, "} }"))
  
  # matern 2, or 3rd order GRF  
  shape = 2
  yacas2 = Ryacas::yacas(paste("prec", shape, ":=MatrixPower(prec0,", shape+1, ")", sep=''))
  yacas2s = Ryacas::yacas(paste("Simplify(prec", shape, ")", sep=''))
  yacasExp = yacas2s$text
  
  yacasChar = gsub("\n", "", as.character(yacasExp))
  
  exprMat = as.matrix(read.table(
          text=gsub("([)][,] |list[(])list[(]", '\n', as.character(yacasChar)), 
          stringsAsFactors=FALSE, sep=','))
  nnIndexMat = as.matrix(NNmat(myGrid, nearest = shape+1))
  precEntriesMat = data.frame(ind=as.vector(nnIndexMat), eqn=as.vector(exprMat))
  precEntriesMat = precEntriesMat[precEntriesMat$ind >0, ]
  precEntriesMat = precEntriesMat[!duplicated(precEntriesMat$ind), ]
  precEntriesMat = precEntriesMat[order(as.integer(precEntriesMat$ind)),]
  precEntriesMat$ind = paste("'", precEntriesMat$ind, "'=", sep='')
  precEntriesMat$last = ','
  precEntriesMat[nrow(precEntriesMat), 'last'] = ''
  
  cat(paste(apply(precEntriesMat, 1, paste, collapse=' '), collapse='\n'), '\n')
  
  # matern 4, or 5th order GRF  
  shape = 4
  yacas4 = Ryacas::yacas(paste("prec", shape, ":=MatrixPower(prec0,", shape+1, ")", sep=''))
  yacas4s = Ryacas::yacas(paste("Simplify(prec", shape, ")", sep=''))
  yacasExp = yacas4s$text
  
  yacasChar = gsub("\n", "", as.character(yacasExp))
  exprMat = as.matrix(read.table(
          text=gsub("([)][,] |list[(])list[(]", '\n', as.character(yacasChar)), 
          stringsAsFactors=FALSE, sep=','))
  nnIndexMat = as.matrix(NNmat(myGrid, nearest = shape+1))
  precEntriesMat = data.frame(ind=as.vector(nnIndexMat), eqn=as.vector(exprMat))
  precEntriesMat = precEntriesMat[precEntriesMat$ind >0, ]
  precEntriesMat = precEntriesMat[!duplicated(precEntriesMat$ind), ]
  precEntriesMat = precEntriesMat[order(as.integer(precEntriesMat$ind)),]
  precEntriesMat$ind = paste("'", precEntriesMat$ind, "'=", sep='')
  precEntriesMat$last = ','
  precEntriesMat[nrow(precEntriesMat), 'last'] = ''
  
  cat(paste(apply(precEntriesMat, 1, paste, collapse=' '), collapse='\n'), '\n')
}