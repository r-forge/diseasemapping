nbMat = rbind(
  c(1,0,0),# self
  c(2,0,1), #up
  c(2,0,-1), #down
  c(2,1,0), #right,
  c(2,-1,0), #left
  c(3,-1,-1),
  c(3,-1,1),
  c(3,1,-1),
  c(3,1,1),
  c(4,0,2), #up
  c(4,0,-2), #down
  c(4,2,0), #right,
  c(4,-2,0), #left
  c(5,-1,-2),
  c(5,-1,2),
  c(5,1,-2),
  c(5,1,2),
  c(5,-2,-1),
  c(5,-2,1),
  c(5,2,-1),
  c(5,2,1),
  c(6,0,3), #up
  c(6,0,-3), #down
  c(6,3,0), #right,
  c(6,-3,0) #left
)


NNmat = function(N,Ny=N, nearest=3, adjustEdges=FALSE) {
  UseMethod("NNmat")	
}

NNmat.Raster = function(N, Ny=N, nearest=3, adjustEdges=FALSE) {
  res = NNmat(ncol(N),nrow(N), nearest, adjustEdges)
  
  attributes(res)$raster= raster(N)
  
  res
}	

NNmat.default = function(N, Ny=N, nearest=3, adjustEdges=FALSE) {
  
  if(nearest<3){
    nbMat = nbMat[nbMat[,1]<=c(2,4)[nearest],] 
  }
  
  Nx = N
  theraster = raster(extent(0,Nx, 0, Ny), nrows = Ny, ncol = Nx)
  
  Ncell = ncell(theraster)
  cellSeq = 1:Ncell
  
  xMat = colFromCell(theraster, cellSeq)
  yMat = rowFromCell(theraster, cellSeq)
  
  xNb = outer(nbMat[,2],xMat, FUN='+')   
  yNb = outer(nbMat[,3],yMat, FUN='+')
  nbCode = matrix(nbMat[,1], nrow(xNb), ncol(xNb))
  
  xNb[union(which(xNb<1),which(xNb>Nx))] = NA
  yNb[union(which(yNb<1),which(yNb>Ny))] = NA
  
  nbIndex = xNb + Nx*(yNb-1)
  nbCol = matrix(1:ncol(xNb), nrow(xNb), ncol(xNb), byrow=TRUE)
  notNa = sort(which(!(is.na(xNb) | is.na(yNb))))
  
  result = sparseMatrix(
    i=nbIndex[notNa],
    j=nbCol[notNa], 
    x=nbCode[notNa]
  )
  
  
  # find id's of cells on the border (for edge correction)
  rasterCells = raster(theraster)
  values(rasterCells) = 1:ncell(theraster)
  
  innerCells = crop(rasterCells, 
    extent(theraster, 
      nearest+1, nrow(theraster)-nearest, 
      nearest+1, ncol(theraster)-nearest))
  
  innerCells = sort(values(innerCells))
  outerCells = sort(setdiff(values(rasterCells), innerCells))
  
  if(adjustEdges) {
    # set entries requiring edge correction to missing
    result[outerCells, outerCells] = NA
  }
  
  result=forceSymmetric(result)
  
  attributes(result)$Nx = Nx
  attributes(result)$Ny = Ny
  attributes(result)$adjustEdges = adjustEdges
  attributes(result)$nearest = nearest
  attributes(result)$raster = theraster
  attributes(result)$cells = list(inner=innerCells, outer=outerCells)
  
  return(result)
}