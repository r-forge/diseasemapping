
x = seq(1, 3, len=2)
(nu = 1.5)
(nuround = round(nu+0.5))
(mu = nu - nuround)

y1 = besselK(x, nu)
y2 = besselK(x, nu, TRUE)
range(log(y2) - x - log(y1), na.rm=T)

x.vec = x
nu.vec = rep(nu, length(x))
jj <- .C("bessel_Knu_scaled_e", as.double(nu.vec), as.double(x.vec), 
    as.integer(length(x.vec)), val = as.double(x.vec), err = as.double(x.vec), 
    status = as.integer(0 * x.vec), PACKAGE = "gsl")
y3 = jj$val

range(y2-y3, na.rm=T)



devtools::load_all("/home/patrick/workspace/diseasemapping/pkg/geostatsgpu")
y4 = .C("bessel_Knu_scaled", nu, as.double(x), as.double(x), 
    as.integer(length(x)), PACKAGE='geostatsgpu')[[3]] 

range(y2-y4, na.rm=TRUE)

y5 = .C("bessel_K_temme", as.integer(nuround), as.double(mu), as.double(x), as.double(x), 
    as.integer(length(x)), PACKAGE='geostatsgpu')[[4]] 

range(y2-y5, na.rm=TRUE)

y6 = .C("bessel_K_p", 
    as.integer(nuround),  as.double(mu),  as.double(nu),
    as.double(x), as.double(x), 
    as.integer(length(x)), PACKAGE='geostatsgpu')[[5]] 
range(y2-y6, na.rm=TRUE)



range(y5-y4, na.rm=TRUE)
range(y5-y6, na.rm=TRUE)
range(besselK(x, nu)*exp(x) - y6)




.C("Rtemme_gamma", 
    as.double(nu), 
    as.double(0), as.double(0), 
    as.double(0), as.double(0), PACKAGE = "geostatsgpu")

