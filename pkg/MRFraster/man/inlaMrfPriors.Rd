\name{inlaMrfPriors}
\alias{inlaMrfPriors}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
inlaMrfPriors(rangePmode = 1, rangeShape = 1, sdPmode = 0.1, varShape = 2, gridsize = 1, cellVar = "cellID")
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{rangePmode}{
%%     ~~Describe \code{rangePmode} here~~
}
  \item{rangeShape}{
%%     ~~Describe \code{rangeShape} here~~
}
  \item{sdPmode}{
%%     ~~Describe \code{sdPmode} here~~
}
  \item{varShape}{
%%     ~~Describe \code{varShape} here~~
}
  \item{gridsize}{
%%     ~~Describe \code{gridsize} here~~
}
  \item{cellVar}{
%%     ~~Describe \code{cellVar} here~~
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
}
\references{
%% ~put references to the literature/web site here ~
}
\author{
%%  ~~who you are~~
}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
##---- Should be DIRECTLY executable !! ----
##-- ==>  Define data, use random,
##--	or do  help(data=index)  for the standard data sets.

## The function is currently defined as
function (rangePmode = 1, rangeShape = 1, sdPmode = 0.1, varShape = 2, 
    gridsize = 1, cellVar = "cellID") 
{
    if (length(grep("^Raster", class(gridsize)))) {
        gridsize = xres(gridsize)
    }
    range = list(prior = "loggamma", param = c(shape = rangeShape, 
        scale = (rangeShape + 1) * rangePmode/gridsize))
    attributes(range)$ci95RangeKm = sort(gridsize/qgamma(c(0.025, 
        0.975), shape = range$param["shape"], scale = 1/range$param["scale"]))
    prec = list(prior = "loggamma", param = c(shape = varShape, 
        scale = sdPmode^2 * (varShape + 1)))
    attributes(prec)$ci95sd = sort(1/sqrt(qgamma(c(0.025, 0.975), 
        shape = prec$param["shape"], scale = 1/prec$param["scale"])))
    return(list(range = range, prec = prec, formula = as.formula(paste(".~.+f(", 
        cellVar, ", model='matern2d',nrow=", attributes(dataForInla)$nrow, 
        ", ncol=", attributes(dataForInla)$ncol, ", nu=2,hyper = list(range=", 
        as.character(thePriors["range"]), ", prec=", as.character(thePriors["prec"]), 
        ") )"))))
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
