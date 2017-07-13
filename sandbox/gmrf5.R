library(Ryacas)
yacas("E4:={ {0,b,0},{b,1,b},{0,b,0} }")

simplify2array(with(list(b=12), eval(yacas("(E4)")$text)))


library('geostatsp')
matrix(NNmat(9, 9, nearest=4)[,9*4+5], 9, 9)

params=c(oneminusar=0.1, shape=0, conditionalVariance=1)

Ngrid = 5
myGrid = squareRaster(extent(-Ngrid,Ngrid,-Ngrid,Ngrid), 1+2*Ngrid)
midcell = cellFromRowCol(myGrid, 
    round(nrow(myGrid)/2), round(ncol(myGrid)/2)) 
cornerSeq = list(row = seq(1,Ngrid+1), col = seq(Ngrid+1, len=Ngrid+1))


# precision matrix without adjusting for edge effects
precMat =maternGmrfPrec(myGrid, param=params) 

attributes(precMat)$info$precisionEntries
attributes(precMat)$info$theo

# the middle cell
precMid=matrix(precMat[,midcell], 
    nrow(myGrid), ncol(myGrid), byrow=TRUE)

precMid[cornerSeq$row, cornerSeq$col]


precMatChar = as.character(precMat)
precMatChar[precMatChar==as.character(1/attributes(precMat)$info$theo['a'])] = '-a'
precMatChar = matrix(precMatChar, nrow(precMat), ncol(precMat))

precForYacas = apply(precMatChar,1,paste,collapse=',')
precForYacas = paste(precForYacas, collapse='},{')
library('Ryacas')
yacas0 =yacas(paste("prec0:= { {", precForYacas, "} }"))

# matern zero
matrix(
          simplify2array(
              with(list(a=1/attributes(precMat)$info$theo['a']), 
                  eval(yacas0$text))
          )[,midcell], 
          nrow(myGrid), ncol(myGrid), byrow=TRUE)[cornerSeq$row, cornerSeq$col]
matrix(
          maternGmrfPrec(
              myGrid, 
              param = c(shape=0, 
                  attributes(precMat)$info$theo[c('conditionalVariance','oneminusar')])
          )[,midcell], 
          nrow(myGrid), ncol(myGrid), byrow=TRUE)[cornerSeq$row, cornerSeq$col]
          
# matern 1

yacas1 = yacas("prec1:=MatrixPower(prec0,2)")

(fromYacas = matrix(
          simplify2array(
              with(list(a=1/attributes(precMat)$info$theo['a']), 
                  eval(yacas1$text))
          )[,midcell], 
          nrow(myGrid), ncol(myGrid), byrow=TRUE))[cornerSeq$row, cornerSeq$col]

(fromGmrfPrec = matrix(
          maternGmrfPrec(
              myGrid, 
              param = c(shape=1, 
                  attributes(precMat)$info$theo[c('conditionalVariance','oneminusar')])
          )[,midcell], 
          nrow(myGrid), ncol(myGrid), byrow=TRUE))[cornerSeq$row, cornerSeq$col]


# matern 2
yacas2 = yacas("MatrixPower(prec0,3)")

(fromYacas = matrix(
          simplify2array(
              with(list(a=1/attributes(precMat)$info$theo['a']), 
                  eval(yacasExp$text))
          )[,midcell], 
          nrow(myGrid), ncol(myGrid), byrow=TRUE))[cornerSeq$row, cornerSeq$col]

(fromGmrfPrec = matrix(
          maternGmrfPrec(
              myGrid, 
              param = c(shape=2, 
                  attributes(precMat)$info$theo[c('conditionalVariance','oneminusar')])
          )[,midcell], 
          nrow(myGrid), ncol(myGrid), byrow=TRUE))[cornerSeq$row, cornerSeq$col]


# matern 4
shape = 4
yacas4 = yacas(paste("MatrixPower(prec0,", shape+1, ")")
yacasExp = yacas4$text

(fromYacas = matrix(
          simplify2array(
              with(list(a=1/attributes(precMat)$info$theo['a']), 
                  eval(yacasExp))
          )[,midcell], 
          nrow(myGrid), ncol(myGrid), byrow=TRUE))

(fromGmrfPrec = matrix(
          maternGmrfPrec(
              myGrid, 
              param = c(shape=shape, 
                  attributes(precMat)$info$theo[c('conditionalVariance','oneminusar')])
          )[,midcell], 
          nrow(myGrid), ncol(myGrid), byrow=TRUE))


# NN 5

 
nnIndices = function(index, baseDist) {
rotDist =  Mod(baseDist)*exp(1i*(Arg(baseDist)+seq(0, 2*pi, len=5)))
rotDist = c(rotDist, Im(rotDist)+1i*Re(rotDist))
rotDist = unique(round(rotDist))

cbind(index, x=Re(rotDist), y=Im(rotDist))
}

rbind(
nnIndices(10,baseDist = 3+2i),
nnIndices(11,baseDist = 4+1i),
nnIndices(12,baseDist = 5+0i)
)
