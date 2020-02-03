#' Matrix of coefficients
#' 
#' @description Produces a matrix of all coeficients suitable for printing
#'
#' @param x model output from lm, glm, lmer, glmer, glmmTMB
#' @param maxChar truncate characters
#' @param link link function
#' @param ... additional arguments for confint
#' @return data frame with estimates and confidence interval.  coefficients are exponentiated if link is log or logit.
#' @export

coefTable = function(x, maxChar = 6, 
  link=family(x)$link, ...) {
#  result = list(confint=confint(x))
result = list(confint = confint(x, ...))
result$tableRaw = summary(x)$coef
if(all(class(x) == 'glmmTMB')) {
  result$table = result$confint
  rownames(result$table ) = gsub("^cond[.]", "", 
    rownames(result$table))
  rownames(result$table ) = gsub("cond[.]Std[.]Dev[.][(]Intercept[)]", "SD", rownames(result$table))

} else {
result$table = cbind(result$tableRaw[,'Estimate', drop=FALSE],
  result$confint[
    match(rownames(result$tableRaw),rownames(result$confint)), ])

otherParameters = setdiff(
  rownames(result$confint), rownames(result$table))

result$table = rbind(result$table, 
  cbind(data.frame(Estimate = rep(NA, length(otherParameters))), 
    result$confint[otherParameters, ,drop=FALSE]))
}

result$table = result$table[, c('Estimate', 
  setdiff(colnames(result$confint), 'Estimate')), drop=FALSE]


variables = attributes(terms(x))$term.labels
variables = intersect(variables, colnames(model.frame(x)))
variablesLevels = lapply(model.frame(x)[variables], 
  function(xx) if(is.logical(xx)){ return(c(TRUE,FALSE))} else {return(levels(xx))})
coefNamePrefix = rep(names(variablesLevels), unlist(lapply(variablesLevels, length)))
coefNames = data.frame(
  variable = coefNamePrefix,
  level = unlist(variablesLevels), stringsAsFactors=FALSE)
rownames(coefNames) = paste0(coefNames$variable, coefNames$level)

if(nrow(coefNames)) {
result$table = cbind(
  as.data.frame(result$table),
  coefNames[match(rownames(result$table), rownames(coefNames)), ])
} else {
  result$table = as.data.frame(result$table, stringsAsFactors=FALSE)
  result$table$variable = rownames(result$table)
  result$table$levels = NA
}

interactions = grep("[:]", rownames(result$table))

for(D in interactions) {
result$table[D, 'level'] = gsub(
paste0('(',paste(unique(coefNames$variable), collapse='|'),')'), 
'', rownames(result$table)[D])
result$table[D, 'variable'] = gsub(
  paste0('(', gsub('[:]', '|', result$table[D,'level']),')'),
    '', rownames(result$table)[D])
}

dontHave = which(is.na(result$table$variable))
outcomeColumns = unique(c('Estimate',colnames(result$confint)))

result$table[dontHave,'variable'] = rownames(result$table)[dontHave]
result$table[dontHave,'level'] = ''
result$table = result$table[,c('variable','level',outcomeColumns)]


interceptLevels = unlist(lapply(variablesLevels, 
  function(xx) c(xx, NA)[1]))
if(any(!is.na(interceptLevels))) {
if(sum(nchar(interceptLevels))>12)
  interceptLevels = trimws(substr(interceptLevels,1,maxChar))
}

if('(Intercept)' %in% rownames(result$table)) {
  result$table['(Intercept)','level'] =
      paste(interceptLevels, collapse=':')
  result$table['(Intercept)','variable'] = 'baseline'
}


if(link[1] %in% c('log','logit')) {
  notSd = grep("[.]SD$|^sigma$", rownames(result$table), invert=TRUE)
  result$table[notSd,outcomeColumns] = 
    exp(result$table[notSd,outcomeColumns])
}
if(link[1] == 'logit') {
  result$table['(Intercept)',outcomeColumns] = 
    result$table['(Intercept)',outcomeColumns] /
    (1+result$table['(Intercept)',outcomeColumns])
  
  result$table['(Intercept)','variable'] = 'ref prob'
}
# standard deviations for lmer stuff
if(any(attributes(class(x))$package == 'lme4')) {
result$var = lme4::VarCorr(x)
result$sd = unlist(lapply(result$var, function(xx) attributes(xx)$stddev))
names(result$sd) = gsub("[.][(]Intercept[)]$", "", names(result$sd))
sigmaRow = grep("^[.]sigma$", rownames(result$table))
if(length(sigmaRow)==1) {
  result$table[sigmaRow,'Estimate'] = attributes(result$var)$sc
  result$table[sigmaRow,c('variable','level')] = c('sd','Obs')
}
for(D in grep("^[.]sig[[:digit:]]+", rownames(result$table), value=TRUE)) {
    Dint = as.integer(gsub("[.]sig", "", D))
    result$table[D,'Estimate'] = 
      attributes(result$var[[Dint]])$stddev[1]
    result$table[match(
      D, rownames(result$table)),c('variable','level')] = 
        c('sd', names(result$var)[Dint])
  }
}

# standard deviations for glmmTMB
if(all(class(x) == 'glmmTMB')) {
  sdRow = grep("[.]SD$|^sigma$", rownames(result$table), value=TRUE)
  result$table[sdRow,'variable'] = 'sd'
  result$table[sdRow,'level'] = gsub("[.]SD$", "", sdRow)
  result$table[sdRow,'level'] = gsub("^sigma$", "obs", 
    result$table[sdRow,'level'])
  ziRow = grep("^zi[.]zi", rownames(result$table), value=TRUE)
  result$table[ziRow,'variable'] = 'zero inf'
  result$table[ziRow,'level'] = gsub(
    "^zi[.]zi[~]", "", ziRow)

if(any(family(x)$family == 'Gamma')) {
  obsRow = grep("^sigma$", rownames(result$table))
  # dispersion parameter, 1/shape, 1/sqrt(shape) = Coef of Var
  result$table[obsRow,outcomeColumns] = 
    sqrt(result$table[obsRow,outcomeColumns])
  result$table[obsRow,'level'] = "sd(Y)/E(Y)"
}

  ziIntercept = grep("^zi[.]zi[~][(]Intercept[)]", 
    rownames(result$table), value=TRUE)
  result$table[ziIntercept,outcomeColumns] = 
    result$table[ziIntercept,outcomeColumns] /
    (1 +   result$table[ziIntercept,outcomeColumns] )
  result$table[ziIntercept,'level'] = 'prob'
}

colnames(result$table) = 
  gsub("Estimate", "est", colnames(result$table))

result
}

