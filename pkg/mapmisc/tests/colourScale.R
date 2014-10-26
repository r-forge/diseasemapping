library('mapmisc')
colourScale(NA,breaks=1:4,style='fixed',labels=c('a','b','c'))

colourScale(NA,breaks=1:4,style='fixed')

colourScale(NULL,breaks=1:4,style='fixed')

colourScale(breaks=1:4,style='fixed')

colourScale(NA,breaks=1:4,style='fixed',labels=c('a','b','c'),col=heat.colors(3), opacity=0.5)


colourScale(x=NA,breaks=1:4,style='unique',labels=c('a','b','c','d'),col=heat.colors(4), opacity=0.5)

colourScale(x=NA,breaks=1:4,style='unique',labels=c('a','b','c','d'),col=t(col2rgb(heat.colors(4))), opacity=0.5)


colourScale(x=NA,breaks=1:4,style='unique',labels=c('a','b','c','d'),col=t(col2rgb(heat.colors(4))), opacity=0.5,exclude=2)

colourScale(x=NA,breaks=1:4,style='unique',labels=c('a','b','c','d'),col=t(col2rgb(heat.colors(4))), opacity=0.5,exclude='a')

colourScale(x=sample(0:5, 12,replace=TRUE),
		breaks=1:4,style='unique',labels=c('a','b','c','d'),
		col=t(col2rgb(heat.colors(4))), opacity=0.5,exclude='a')


colourScale(x=sample(0:5, 12,replace=TRUE),
		breaks=1:4,style='unique',labels=c('missing','a','b','c','d','e'),
		col=t(col2rgb(heat.colors(4))), opacity=0.5,exclude='a')


colourScale(seq(0,10),breaks=4,style='equal')
colourScale(seq(0,10),breaks=4,style='quantile')
colourScale(x=seq(0,10),breaks=4,style='equal',exclude=0)
colourScale(x=seq(0,10),breaks=4,style='equal',exclude=c(0,10))
colourScale(x=seq(0,10),breaks=4,style='equal',exclude='nothing')

myraster = raster(matrix(0:8, 3, 3))
colourScale(x=myraster,breaks=4,style='equal')
colourScale(x=myraster,breaks=4,style='quantile')
colourScale(x=myraster,breaks=4,style='unique')
colourScale(x=myraster,breaks=1:4,style='unique',labels=c('a','b','c','d'))
colourScale(x=myraster,breaks=1:4,style='fixed')

colourScale(x=myraster,breaks=4,style='equal',exclude=0)
colourScale(x=myraster,breaks=4,style='unique',exclude=0)

colourScale(x=factor(sample(0:4,10,replace=TRUE)),breaks=4,style='thisShouldBeIgnored',exclude='0')

