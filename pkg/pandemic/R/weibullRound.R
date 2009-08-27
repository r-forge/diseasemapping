censorweibull=function(a,b,cen, params)
{
weibull=0
for(i in 1:length(cen))
{
weibull[i]=0
while(weibull[i]<=cen[i])
{
weibull[i]=round(rweibull(1,a,b))
}
}
weibull
}



rweibullRound = function(N, params=c(shape=1, scale=1)) {

  round(rweibull(N, shape=params["shape"], scale=params["scale"]))
}


dweibullRound = function(x, params=c(shape=1, scale=1)) {

  if(length(x)==0) return(NULL)

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
  prob[x==0]=params["zeros"]+prob [x==0]                                                                 
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

sumWeibullzero=function(k,params1,params2)
{
x=0
for(j in 0:k)
{
for(i in 0:(k-j))
{
x=x+dweibullzero(j,params1)*dweibullzero(i,params2)
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