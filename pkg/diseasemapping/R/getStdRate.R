getStdRate= function(relativeRate, model, referencePopulation, scale=100000) {

relativeRate * sum(predict(model, referencePopulation)) * scale

}