
pcPriorRange = function(q, p=0.5, cellSize=1) {
  if(!is.numeric(cellSize)) {
    cellSize = xres(cellSize)
  } 
  # suppose lower 5% quantile of rangeKm is 50km
# ... of range is 50000/cellSize
# upper 5% quantile (lower 95%) of scale is  cellSize / 50000 
  # 1-exp(-lambda x) = 1-p
  # 1-0.05 = 1-exp(-lambda * cellSize/50000)
# log(0.05) = -lambda * cellSize / 50000
# lambda = -log(0.05) * 50000/cellSize
  # lambda is a multiplicative parameter
  # I call lambda a scale, wikipedia and R call it rate
  quantileCells = q/cellSize
  lambda = -log(p) * quantileCells / log(2)
  resultInla = paste0(
    "expression:
    lambda = ", lambda, ";
    scale = exp(-log_range);
    logdens = -lambda*scale + log(lambda);
    log_jacobian = - log_range;
    return(logdens + log_jacobian);", sep='')

  result= list(
    string = paste0("prior=\"", resultInla, "\"" ),
    param = c(q = q, prob = p), 
    extra = list(
      lambda=unname(lambda),
      medianRange = unname(cellSize*lambda), 
      medianScale = unname(1/(cellSize*lambda)),
      meanScale = unname(1/(log(2)*cellSize*lambda))),
    info = 'exponential prior for scale (pc prior)',
    cellSize = unname(cellSize))

    result$dprior = list(
      range = eval(parse(text=paste0(
        'function(x) x^(-2)*stats::dexp(1/x, rate=', 
        1/result$extra$meanScale, ')'))),
      scale = eval(parse(text=paste0(
        'function(x) stats::dexp(x, rate=', 1/result$extra$meanScale, ')')))
      )
    environment(result$dprior$range) = emptyenv()
    environment(result$dprior$scale) = emptyenv()
  result
}