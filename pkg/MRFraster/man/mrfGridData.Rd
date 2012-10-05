\name{mrfGridData}
\alias{mrfGridData}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
mrfGridData(points, Napprox = 11, covariates = NULL, offset = NULL, buffer = 2)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{points}{
%%     ~~Describe \code{points} here~~
}
  \item{Napprox}{
%%     ~~Describe \code{Napprox} here~~
}
  \item{covariates}{
%%     ~~Describe \code{covariates} here~~
}
  \item{offset}{
%%     ~~Describe \code{offset} here~~
}
  \item{buffer}{
%%     ~~Describe \code{buffer} here~~
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
function (points, Napprox = 11, covariates = NULL, offset = NULL, 
    buffer = 2) 
{
    thebb = bbox(points)
    Nx = Napprox
    cellSize = (thebb[1, 2] - thebb[1, 1])/Nx
    Ny = (thebb[2, 2] - thebb[2, 1])/cellSize
    thebb[2, 2] = thebb[2, 1] + cellSize * Ny
    rasterTemplate = raster(nrows = Ny, ncols = Nx, xmn = thebb[1, 
        1], xmx = thebb[1, 2], ymn = thebb[2, 1], ymx = thebb[2, 
        2])
    count = rasterize(thesim$points, rasterTemplate, field = rep(1, 
        length(thesim$points)), fun = "sum")
    count[][is.na(count[])] = 0
    result = data.frame(count = count[], offset = 2 * log(cellSize))
    if (!is.null(offset)) {
        result$offset = result$offset + resample(offset, rasterTemplate)[]
    }
    if (!is.null(covariates)) {
        covariates@layernames[covariates@layernames == ""] = paste("cov", 
            seq(1, sum(covariates@layernames == "")), sep = "")
        for (Dcov in covariates@layernames) {
            result[[Dcov]] = resample(covariates[[Dcov]], rasterTemplate)[]
        }
    }
    Nx2 = Nx + 2 * buffer
    Ny2 = Ny + 2 * buffer
    xvec = seq(from = buffer + 1, by = 1, len = Nx)
    xvec = rep(xvec, Ny)
    yvec = seq(from = buffer + 1, by = 1, len = Ny)
    yvec = rep(yvec, rep(Ny, Nx))
    index = xvec + Nx2 * (yvec - 1)
    result$cellID = index
    rasterTemplateBig = raster(nrows = Ny2, ncols = Nx2, xmn = thebb[1, 
        1] - buffer * cellSize, xmx = thebb[1, 2] + buffer * 
        cellSize, ymn = thebb[2, 1] - buffer * cellSize, ymx = thebb[2, 
        2] + buffer * cellSize, crs = CRS(proj4string(points)))
    attributes(result)$rasterTemplate = rasterTemplateBig
    attributes(result)$buffer = buffer
    attributes(result)$bbox = thebb
    attributes(result)$nrowForInla = ncol(rasterTemplateBig)
    attributes(result)$ncolForInla = nrow(rasterTemplateBig)
    result
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
