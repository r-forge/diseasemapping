\name{SurvToInla}
\alias{SurvToInla}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Convert Surv objects to survival responses for INLA}
\description{
INLA requires a data frame for survival models
}
\usage{
SurvToInla(x)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{A Surv object}
}
\value{
A data frame with the following columns
  \item{event}{0 (right censored), 1 (observed event), 2 (left censored), or 3 (interval censored)}
  \item{time}{When event=1, this is the time of the event}
  \item{truncation}{The likelihood should condition on the event time being after the truncation time}
  \item{lower}{For right censored or interval censored data, the lower bound for the event time}
  \item{upper}{For left censored or interval censored data, the upper bound for the event time}
}
\author{Patrick Brown}
\seealso{Surv}
\examples{

# right censored
x = Surv(runif(4, 0, 10),  c(0,1,1,0) )
x
SurvToInla(x)

# left censored
x = Surv(runif(4, 0, 10),  c(0,1,1,0), type="left" )
x
SurvToInla(x)


# interval censored data
x = Surv(runif(4, 0, 5), runif(4,5,10),   c(0,1,2,3), type="interval" )
x
SurvToInla(x)

# left truncated, right censored data
x = Surv(runif(4, 0, 5), runif(4,5,10),   c(0,1,1,0))
x
SurvToInla(x)

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
