getInits = function() { 

scale = 1.5
SDscale = 4

result = list()

result[["intercept"]] = sign(startingValues[["intercept" ]]) *
    runif(length(startingValues[["intercept" ]]),
       abs(startingValues[["intercept"]])/scale,
       scale * abs(startingValues[["intercept"]]))


result[["SDCSDUID"]] = sqrt(runif(1,
       startingValues$vars[["CSDUID"]]/scale,
       startingValues$vars[["CSDUID"]]*scale))

result[["RCSDUID"]] = rnorm(length(startingValues[["RCSDUID"]]),
        startingValues[["RCSDUID"]], startingValues$vars[["CSDUID"]]/SDscale)

result[["TCSDUIDSpatial"]] = 1/(sqrt(runif(1,
       startingValues$vars[["CSDUIDSpatial"]]/scale,
       startingValues$vars[["CSDUIDSpatial"]]*scale)))^2

result[["RCSDUIDSpatial"]] = rnorm(length(startingValues[["RCSDUIDSpatial"]]),
        startingValues[["RCSDUIDSpatial"]], startingValues$vars[["CSDUIDSpatial"]]/SDscale)


return(result)

}