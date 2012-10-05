\name{inla.collect.mcmc}
\alias{inla.collect.mcmc}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
inla.collect.mcmc(x, all = T, nBurn = 100)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{
%%     ~~Describe \code{x} here~~
}
  \item{all}{
%%     ~~Describe \code{all} here~~
}
  \item{nBurn}{
%%     ~~Describe \code{nBurn} here~~
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
function (x, all = T, nBurn = 100) 
{
    covs = colnames(x$model.matrix)[-1]
    wd.mcmc = x$.args$working.directory
    dir.mcmc = paste(wd.mcmc, "/", "results.files", sep = "")
    dd = dir(dir.mcmc)
    if (all) {
        result = summarys = as.list(dd)
        names(result) = names(summarys) = dd
        for (id in 1:length(dd)) {
            d.mcmc = paste(dir.mcmc, "/", dd[id], sep = "")
            fnm = paste(d.mcmc, "/", "trace.dat", sep = "")
            if (file.exists(fnm)) {
                allS = read.table(fnm)
                result[[id]] = data.frame(allS[-c(1:nBurn), ])
            }
            if (class(result[[id]]) == "data.frame") {
                tempData = data.frame(result[[id]][-(1:nBurn), 
                  ])
                temp = data.frame(t(rbind(rbind(apply(result[[id]], 
                  2, summary), apply(result[[id]], 2, quantile, 
                  c(0.05, 0.95))), sd = apply(result[[id]], 2, 
                  sd))))
                rownames(temp) = NULL
                summarys[[id]] = temp[, c(4, 9, 7, 8, 1, 2, 3, 
                  5, 6)]
            }
        }
    }
    lastConf = as.list(1:3)
    names(lastConf) = c("idTable", "Theta", "x")
    lastDIR = paste(dir.mcmc, "/", dd[grep("^last", dd, ignore.case = T)], 
        sep = "")
    fnmTable = paste(lastDIR, "/", "idx-table.dat", sep = "")
    fnmTheta = paste(lastDIR, "/", "theta.dat", sep = "")
    fnmX = paste(lastDIR, "/", "x.dat", sep = "")
    if (file.exists(fnmTable)) {
        lastConf[["idTable"]] = read.table(fnmTable)
    }
    if (file.exists(fnmTheta)) {
        lastConf[["Theta"]] = read.table(fnmTheta)
    }
    if (file.exists(fnmX)) {
        lastConf[["x"]] = read.table(fnmX)
    }
    offDir = paste(dir.mcmc, "/", dd[grep("^totaloffset", dd, 
        ignore.case = T)], sep = "")
    offTable = paste(offDir, "/", "totaloffset.dat", sep = "")
    if (file.exists(offTable)) {
        totalOffset = read.table(offTable)
    }
    if (all) {
        x = list(result, summarys, lastConf, totalOffset)
        names(x) = c("samples", "summary", "lastConf", "offset")
        torename = paste("^fixed\\.effect(0+)", seq(from = 2, 
            length = length(covs)), "$", sep = "")
        names(torename) = covs
        torenamehyp = paste("^hyperparameter([[:digit:]]|-|random\\.effect)+parameter", 
            rep(0:1, c(2, 2)), rep(c("", "-user-scale"), 2), 
            "$", sep = "")
        names(torenamehyp) = paste(rep(c("log", ""), 2), rep(c("Precision", 
            "Range"), c(2, 2)), sep = "")
        torename = c(torename, intercept = "^fixed\\.effect(0+)1$", 
            torenamehyp)
        for (D1 in c("samples", "summary")) {
            for (D2 in 1:length(torename)) names(x[[D1]]) = gsub(torename[D2], 
                names(torename)[D2], names(x[[D1]]))
        }
    }
    x
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
