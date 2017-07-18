#' Matrix for confidence intervals
#' 
#' @description Produces a matrix suitable for multiplying by results of summary or predict functions to give confidence intervals of a desired quantile
#'
#' @param p coverage of confidence interval
#' @param se.fit row names of result are 'fit' and 'se.fit'
#' @return matrix with three columns (estimate, upper and lower bound of CI) and two rows
#' @examples
#' 
#' (myCiMat = ciMat(0.8))
#' 
#' clotting <- data.frame(
#'     u = c(5,10,15,20,30,40,60,80,100),
#'     lot1 = c(118,58,42,35,27,25,21,19,18),
#'     lot2 = c(69,35,26,21,18,16,13,12,12))
#' glmRes = stats::glm(lot1 ~ log(u), data = clotting, family = Gamma)
#' # CI on the natural scale
#' exp(summary(glmRes)$coef[,rownames(myCiMat)] %*% myCiMat)
#' 
#' (myCiMatPred = ciMat(0.99, se.fit=TRUE))
#' glmPred = do.call(cbind,stats::predict(glmRes, se.fit=TRUE))
#' exp(glmPred[,rownames(myCiMatPred)] %*% myCiMatPred)
#' 
#' @export
ciMat = function(p=0.95, se.fit = FALSE) {
  p = (1-p)/2
  p = sort(unique(c(p, 1-p)))
  res = rbind(
      1,
      c(0, qnorm(p))
  )
  colnames(res) = c('est', format(100*p))
  if(se.fit) {
    rownames(res) = c('fit','se.fit')
  } else {
    rownames(res) = c('Estimate','Std. Error')
  }
  res
}