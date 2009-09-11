\name{expandCensus}
\Rdversion{1.1}
\alias{expandCensus}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
expandCensus(popdata, year.range = range(as.integer(names(popdata))), years = as.integer(names(popdata)))
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{popdata}{
%%     ~~Describe \code{popdata} here~~
}
  \item{year.range}{
%%     ~~Describe \code{year.range} here~~
}
  \item{years}{
%%     ~~Describe \code{years} here~~
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
function(popdata,year.range=range(as.integer(names(popdata))),years=as.integer(names(popdata))){


  times <- c(year.range[1], sort(years), year.range[2])
  times <- as.numeric(times)
  inter <- diff(times)/2
  nseq <- 1:length(inter) - 1
  mseq <- 2:length(inter)
  interval <- inter[mseq] + inter[nseq]
  interval[1] <-interval[1] + inter[1]
  interval[length(interval)] <-interval[length(interval)] + inter[length(inter)]
  interval<-ceiling(interval)

  newp<-list()
  index=1
  for (k in interval){
  
    if(index==1) n<-1 else n <- length(newp)+1
    for (i in n:(n+k-1)){
      newp[[i]]<-popdata[[index]]
    }
  index=index+1
 }
 
  names(newp) = year.range[1]:year.range[2]


newp
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
