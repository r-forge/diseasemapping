
x = raster(nrows=11,ncols=11,xmn=0,xmx=10,ymn=0,ymx=10)

temp=GaussRF(x, model="whittle", 
		param=c(mean=0, variance=1, nugget=0, 
				scale=2, alpha=2))


x = SpatialPoints(cbind(1:4, 11:14))

 temp=GaussRF(x, model="whittle", 
		param=c(mean=0, variance=1, nugget=0, 
				scale=2, alpha=2))


x = cbind(1:4, 1:4)
temp=GaussRF(x, model="whittle", 
		param=c(mean=0, variance=1, nugget=0, 
				scale=2, alpha=2))
