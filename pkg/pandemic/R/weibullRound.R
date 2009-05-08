rweibullRound = function(N, params=c(shape=1, scale=1)) {

  round(rweibull(N, shape=params["shape"], scale=params["scale"]))
}


dweibullRound = function(x, params=c(shape=1, scale=1)) {
  xLow = pmax(x-0.5, 0)
  xHigh = x+0.5
  pweibull(xHigh, shape=params["shape"], scale=params["scale"]) - 
    pweibull(xLow, shape=params["shape"], scale=params["scale"])
                                                                   

}

pweibullRound=function(x, params) {
  pweibull(x+0.5, shape=params["shape"], scale=params["scale"])
}