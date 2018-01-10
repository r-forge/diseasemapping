
breaksForRates = function(
    x,
    breaks = 10,
    transform = 0.1,
    multiples = c(2, 4, 5, 10)) {
  
  if(methods::existsMethod(raster::maxValue, class(x)))
    x = raster::maxValue(x)
  theMax = max(x)
  
  if(is.character(transform))
    transform = c(sqrt=0.5, none=1)[transform]
  
  SpropOfMax = multiples
  
  DpropOfMax = 0
  theBreaksP1 = NULL
  maxIter = 50
  Diter = 0
  toRound = Inf

  while(Diter < maxIter) {
    toRoundOld = toRound
    Diter = Diter + 1
    DpropOfMax = DpropOfMax + 1
    theBreaks = theBreaksP1
    if(DpropOfMax > length(SpropOfMax)) {
      SpropOfMax = SpropOfMax*10
      DpropOfMax = 1 
    }
    roundPropOfMax = SpropOfMax[DpropOfMax]
    
    toRound =   10^round(log10(theMax))/roundPropOfMax
    breakSeq1 = seq(0, theMax^(transform), len=breaks)^(1/transform)
    breakSeq1[length(breakSeq1)] = theMax +.Machine$double.eps
    ceilingSeq = ceiling(breakSeq1/toRound)
    
    moreRounding = 10^floor(log10(ceilingSeq)-1)
    
    breakSeq = ceiling(ceilingSeq/moreRounding)*moreRounding*toRound
    breakSeq[1] = 0
    breakSeq = sort(unique(signif(breakSeq,10)))

    theBreaksP1 = sort(unique(c(
      breakSeq[breakSeq < 2*toRoundOld], 
      theBreaksP1)))


#    theBreaksP1 = sort(unique(signif(breakSeq,10)))
    
    if(length(theBreaksP1) >= breaks) {
      break()
    }
    
  }
  
  if(length(theBreaksP1) > (2+breaks) ) 
    theBreaksP1 = theBreaks
  
  if(sum(theBreaksP1 > theMax) > 2) {
    theBreaksP1 = theBreaksP1[-length(theBreaksP1)]
  }

  theBreaksP1
}