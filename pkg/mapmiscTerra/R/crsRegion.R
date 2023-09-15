
crsRegionEllipse = function(x, offset) {
  # need this here because of unknown global variable warning in package checks
  isoSp = vect(mapmisc::isohedron, crs=crsLL)
  isoT = project(isoSp, x)
  if(missing(offset)) {
    offset = apply(matrix(ext(isoT), ncol=2), 2, mean)  
  }
  isoT = crds(isoT) 
  isoT = isoT - matrix(offset, nrow(isoT), 2, byrow=TRUE)


  isoTc = isoT[,1] + 1i*isoT[,2]
  isoTc = cbind(Arg(isoTc), Mod(isoTc))
#  isoTc = isoTc[isoTc[,2] > quantile(isoTc[,2], 0.1),]

  axisX = max(isoTc[,2])

  rfun = function(theta, a, b) sqrt((b*cos(theta))^2 + (a*sin(theta))^2)
  objFun = function(param, x=isoTc) {
    100*sum(pmax(0, x[,2] -rfun(x[,1], param[1], param[2]))^2) + sum(log(param))
  }
  oo = stats::optim(c(axisX, axisX), objFun, x=isoTc)


  angleSeq = seq(0, 2*pi, len=1001)
  ellipseP = cbind(cos(angleSeq)*oo$par[2], sin(angleSeq)*oo$par[1])

  ellipseP + matrix(offset, nrow(ellipseP), 2, byrow=TRUE)
}