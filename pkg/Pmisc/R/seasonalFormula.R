#' Seasonal formula
#' 
#' @description A formula with sines and cosines
#'
#' @param period 1/frequency, 365.25 for daily data with yearly cycles
#' @param harmonics multiplied by period, one sin/cosine pair per harmonic
#' @param var variable name
#' @return a formula with sine/cosine pairs
#' @example seasonalFormula(period=365.25, harmonics=1:2, var='x')
#' @example ff1 <- seasonalFormula(period=12, harmonics=1:2, var='x')
#' @example ff2 <- update.formula(y~x, ff1)
#' @example model.matrix(ff2, data.frame(y=1:4,x=1:4))
#' @export

seasonalFormula=function(period=365.25,harmonics=1:2,var='x') {

	periodOverHarmonics = format(period/(2*harmonics), digits=10)

	res1 = paste0('(', var, '/', periodOverHarmonics, ')')
	res2 = paste(outer(c('sinpi','cospi'), res1, FUN=paste0), collapse=' + ')

	stats::as.formula(paste('.~.+', res2), env=new.env())
}

