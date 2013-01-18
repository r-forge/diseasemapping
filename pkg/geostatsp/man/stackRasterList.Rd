\name{stackRasterList}
\alias{stackRasterList}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
stackRasterList(x, template = x[[1]], method = "ngb")
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{
%%     ~~Describe \code{x} here~~
}
  \item{template}{
%%     ~~Describe \code{template} here~~
}
  \item{method}{
%%     ~~Describe \code{method} here~~
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
function (x, template = x[[1]], method = "ngb") 
{
    if (length(method) == length(x)) {
        if (all(names(x) \%in\% names(method))) 
            method = method[names(x)]
    }
    else {
        method = rep(method, length(x))
    }
    names(method) = names(x)
    result = template
    template = as(template, "BasicRaster")
    for (D in names(x)) {
        if (as(x[[D]], "BasicRaster") == template) {
            result = addLayer(result, x[[D]])
        }
        else {
            result = addLayer(result, projectRaster(x[[D]], template, 
                method = method[[D]]))
        }
    }
    if (length(x) == (length(names(result)) + 1)) 
        result = result[[-1]]
    names(result) = names(x)
    result
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
