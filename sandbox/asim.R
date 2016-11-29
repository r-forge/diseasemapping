

## Equation 23 ------------------------------------------------------------------------
library(rpanel)

my.panel <- function(panel){
  ## parameters
  lambda <- panel$lambda
  sigma2 <- panel$sigma2
  phi2 <- panel$phi2
  mu <- panel$mu
  
  ## Equation 23
  cuv <- function(dist, lambda, phi2, sigma2, mu){
    u <- dist[1]
    v <- dist[2]
    exp(-phi2*(u^2)/2)/((lambda + (u^2)*sigma2/2)^2 + (v + u*mu)^2)
  }
  
  ## build perspective vector
  space <- seq(-5, 5, length=50)
  time <- seq(-5, 5, length=50)
  spacetime <- expand.grid("space"=space, "time"=time)
  
  z <- matrix(apply(as.matrix(spacetime), 1, cuv,
          lambda = panel$lambda, sigma2 = panel$sigma2,
          phi2 = panel$phi2, mu = panel$mu
  		),
  		length(space), length(time))
  
  ## Contour plot
  contour(x = space, y = time, z = z,
      xlab = "Space, u",
      ylab = "Time, v",
      main = "Surface of C(u, v)")
  
  ## If you want to see the perspective, 
  ## please uncomment the following lines and 
  ## comment the contour lines above
  # persp(x = space,
  #       y = time,
  #       z = z,
  #       phi = 20,
  #       theta = 30,
  #       xlab = "Space, u",
  #       ylab = "Time, v",
  #       zlab = "C(u, v)",
  #       main = "Surface of C(u, v)"
  # )
  panel
}


x11()
panel <- rp.control(interval=c(0.1, 3))

## Equation 23
rp.slider(panel, lambda,
    from=0.001, to=10, initval=1,
    showvalue=TRUE,
    action=my.panel,
    #resolution = 0.1,
    title="lambda")

rp.slider(panel, sigma2,
    from=0.001, to=10, initval=1,
    showvalue=TRUE,
    action=my.panel,
    #resolution = 0.1,
    title="sigma2")

rp.slider(panel, mu,
    from=-1.5, to=1.5, initval=0,
    showvalue=TRUE,
    action=my.panel,
    #resolution = 0.1,
    title="mu")

rp.slider(panel, phi2,
    from=0.001, to=10, initval=1,
    showvalue=TRUE,
    action=my.panel,
    #resolution = 0.1,
    title="phi2")


## Equation 35 ------------------------------------------------------------------------
library(rpanel)

my.panel <- function(panel){
  ## parameters
  lambda <- panel$lambda
  sigma <- panel$sigma
  phi <- panel$phi
  mu <- panel$mu
  a <- panel$a
  b <- panel$b
  
  ## Equation 35
  chk <- function(dist, lambda, sigma, phi, mu, a, b){
    h <- dist[1]
    k <- dist[2]
    gamma <- sqrt(2*lambda)
    g1 <- pnorm((sigma*mu*k - gamma*b^2 - sigma*h)/(sigma*b), 0, 1)
    g2 <- pnorm((sigma*h - gamma*b^2 - sigma*mu*k)/(sigma*b), 0, 1)
    e1 <- exp(((gamma^2) * (phi^2))/(2*sigma^2))
    (1/(2*gamma*sigma))*e1*(exp(a)*g1 + exp(-a)*g2)
  }
  
  ## build perspective vector
  space <- seq(-5, 5, length=50)
  time <- seq(-5, 5, length=50)
  spacetime <- expand.grid("space"=space, "time"=time)
  
  z <- matrix(apply(as.matrix(spacetime), 1, chk,
          phi = panel$phi, lambda = panel$lambda, sigma = panel$sigma,
          mu = panel$mu, a = panel$a, b = panel$b
      ),
  		length(space), length(time))
  
  ## plot
  contour(x = space, y = time, z = z,
      xlab = "Space, h",
      ylab = "Time, k",
      main = "Surface of C(h, k)")
  
  ## If you want to see the perspective, 
  ## please uncomment the following lines and 
  ## comment the contour lines above
  # persp(x = space,
  #       y = time,
  #       z = z,
  #       phi = 20,
  #       theta = 30,
  #       xlab = "Space, h",
  #       ylab = "Time, k",
  #       zlab = "C(h, k)",
  #       main = "Surface of C(h, k)"
  # )
  panel
}


x11()
panel <- rp.control(interval=c(0.1, 3))

## Equation 35
rp.slider(panel, lambda,
    from=0.001, to=10, initval=1,
    showvalue=TRUE,
    action=my.panel,
    #resolution = 0.1,
    title="lambda")

rp.slider(panel, sigma,
    from=0.001, to=10, initval=1,
    showvalue=TRUE,
    action=my.panel,
    #resolution = 0.1,
    title="sigma")

rp.slider(panel, mu,
    from=-1.5, to=1.5, initval=0,
    showvalue=TRUE,
    action=my.panel,
    #resolution = 0.1,
    title="mu")

rp.slider(panel, phi,
    from=0.001, to=10, initval=1,
    showvalue=TRUE,
    action=my.panel,
    #resolution = 0.1,
    title="phi")

rp.slider(panel, a,
    from=-1.5, to=1.5, initval=0.5,
    showvalue=TRUE,
    action=my.panel,
    #resolution = 0.1,
    title="a")

rp.slider(panel, b,
    from=0.001, to=10, initval=1,
    showvalue=TRUE,
    action=my.panel,
    #resolution = 0.1,
    title="b")