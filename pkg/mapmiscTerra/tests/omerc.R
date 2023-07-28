
library('mapmisc')
data('netherlands')

nldCities = unwrap(nldCities)


x=nldCities
xbox = as.polygons(ext(nldCities), crs = crs(nldCities))

bob=function(angle, ...){
  y = project(x, omerc(x, angle, ellipse=FALSE, ...))
  nld2 = project(xbox, crs(y))
  map.new(nld2)
  abline(v=0, col='grey')
  abline(h=0, col='grey')
  plot(nld2,add=TRUE, lwd=2, lty=3)
  plot(ext(y), add=TRUE, col='orange')
  plot(y,cex=0.2, col='red',add=TRUE)
  text(y,labels=y$name, cex=0.5, col=col2html('blue',0.4))
   mtext(paste(angle,collapse=' '),side=3,outer=FALSE,line=-1)
   scaleBar(y,'topright')
   return(invisible(crs(y)))
  
 }



if(!interactive()) pdf("omerc.pdf")

  par(mfrow=c(3,3))


  	bob(0)
  
  bob(89)
  
  bob(45)
  
  bob(-45)

  bob(180-45)

  
  bob(-180-45, post=-45)
  
  bob(45, post='north')
  bob(45, post='wide')
  bob(45, post='tall')
  
  


par(mfrow=c(3,2))
  bob((-10):10)
  
  bob(seq(-170,-190))
  
  bob((-10):10, post='north')
  bob((-10):10, post='wide')
  bob((-10):10, post='long')
  
  bob((-10):10, post=90)

N = 12
somePoints = vect(
    cbind(runif(N,-5,40), runif(N,40,70)),
    atts=data.frame(name=1:N),
    crs=crsLL
    )
    x=somePoints
    xbox = as.polygons(ext(somePoints), crs(somePoints))
 
  par(mfrow=c(3,2))
  
  
  bob((-10):10, preserve=x, post=10)
  bob((-10):10, preserve=x, post='none')
  
  bob(seq(-170,-190), preserve=x)
  
  scaleBar(bob((-10):10, preserve=x, post='north'), c(0,0), seg.len=0)

  
  scaleBar(bob((-10):10, preserve=x, post='tall'), c(0,0), seg.len=0)
  
  scaleBar(bob((-10):10, preserve=x, post='wide'), c(0,0), seg.len=0)
  

if(!interactive())   dev.off()

