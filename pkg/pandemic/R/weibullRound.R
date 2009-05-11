rweibullRound = function(N, params=c(shape=1, scale=1)) {

  round(rweibull(N, shape=params["shape"], scale=params["scale"]))
}


dweibullRound = function(x, params=c(shape=1, scale=1)) {
  xLow = pmax(x-0.5, 0)
  xHigh = x+0.5
  pweibull(xHigh, shape=params["shape"], scale=params["scale"]) - 
    pweibull(xLow, shape=params["shape"], scale=params["scale"])
                                                                   

}



dweibullzero= function(x, params=c(shape=1, scale=1,zeros=0)) {
  xLow = pmax(x-0.5, 0)
  xHigh = x+0.5
 prob= pweibull(xHigh, shape=params["shape"], scale=params["scale"]) - 
    pweibull(xLow, shape=params["shape"], scale=params["scale"])
 prob=(1-params["zeros"])*prob
  if(x==0) prob=params["zeros"]+prob                                                                  
 prob
}

pweibullRound=function(x, params) {
  pweibull(x+0.5, shape=params["shape"], scale=params["scale"])
}

sumWeibull=function(k,params1,params2)
{
x=0
for(j in 0:k)
{
for(i in 0:(k-j))
{
x=x+dweibullRound(j,params1)*dweibullRound(i,params2)
}
}
x
}

simWeibullzero=function(params)
{
u=runif(1)
x=-1
while(u>0)
{
x=x+1
u=u-dweibullzero(x,params)
}
x
}