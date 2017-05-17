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
  cellSeq = values(theraster) = 1:Ncell

  # find id's of cells on the border (for edge correction)  
  innerCells = crop(theraster, 
    extend(extent(theraster), -nearest))
  
  innerCells = sort(values(innerCells))
  outerCells = sort(setdiff(values(theraster), innerCells))
  
  
  xMat = colFromCell(theraster, cellSeq)
  yMat = rowFromCell(theraster, cellSeq)

  xNb = outer(nbMat[,2],xMat, FUN='+')   
  yNb = outer(nbMat[,3],yMat, FUN='+')

  # which neighbour pairs are outside of the grid
  nbOutsideRegion = which(Matrix(xNb < 1) | Matrix(xNb > Nx) | 
    Matrix(yNb < 1) | Matrix(yNb> Ny))

  xNb = xNb[-nbOutsideRegion]
  yNb = yNb[-nbOutsideRegion]
  
  nbPairs = cbind(
  i = xNb + Nx*(yNb-1),
  j = matrix(cellSeq, nrow(nbMat), Ncell, byrow=TRUE)[-nbOutsideRegion],
  x = matrix(nbMat[,1], nrow(nbMat), Ncell)[-nbOutsideRegion])
  
  if(adjustEdges) {
    # set entries requiring edge correction to missing
    # allocating the memory here instead of in maternGmrfPrec is faster
    outerCellsPairs = expand.grid(i=outerCells, j=outerCells, x=NA)
    nbPairs = rbind(
      nbPairs, outerCellsPairs
      )
  }

  nbPairs = as.list(as.data.frame(nbPairs))
  nbPairs$dims = c(Ncell, Ncell)
  nbPairs$symmetric = TRUE
  nbPairs$dimnames = list(
    cellSeq, cellSeq
    )
  nbPairs$check = FALSE
    
  result = do.call(sparseMatrix, nbPairs)
  
  attributes(result)$Nx = Nx
  attributes(result)$Ny = Ny
  attributes(result)$adjustEdges = adjustEdges
  attributes(result)$nearest = nearest
  attributes(result)$raster = theraster
  attributes(result)$cells = list(inner=innerCells, outer=outerCells)
  
  return(result)
}