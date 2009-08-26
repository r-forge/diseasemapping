\name{pandemicParams}
\Rdversion{1.1}
\alias{pandemicParams}
\title{
Paramter Generators
}
\description{
\code{pandemicParams} is a generic function for generating paramters for simulations.
}
\usage{
pandemicParams(InfOns = c(mean = 1, shape = 1, zeros = 0.1), OnsMedM = c(mean = 3, shape = 2, zeros = 0.1), OnsMedS = c(mean = 2.5, shape = 1, zeros = 0.1), OnsMedD = c(mean = 2, shape = 1, zeros = 0.1), MedRec = c(mean = 7, shape = 2, zeros = 0.1, lost = 0.9), MedHospS = c(mean = 2.5, shape = 2, zeros = 0.1), MedHospD = c(mean = 1.5, shape = 1.5, zeros = 0.1), HospRec = c(mean = 4, shape = 1.5, zeros = 0.1), HospDeath = c(mean = 3, shape = 1, zeros = 0.1), probs = c(M = 0.6, S = 0.3, D = 0.1),
 ageProbs=NULL)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{InfOns}{time from stage infection to onsets.}
  \item{OnsMedM}{time from stage onsets to on medication for mild cases.}
  \item{OnsMedS}{time from onsets to on medication for serious cases.}
  \item{OnsMedD}{time from onsets to on medication for deadly cases.}
  \item{MedRec}{time from on medication to recover.}
  \item{MedHospS}{time from on medication to hospitalization for serious cases.}
  \item{MedHospD}{time from on medication to hospitalization for deadly cases.}
  \item{HospRec}{time from hospitalization to recover.}
  \item{HospDeath}{time from hospitalization to death.}
  \item{probs}{probabily of mild cases, serious cases and deadly cases.}
\item{ageProbs}{probabilities which vary with age. A list with elements S and/or D. A probability can be included in only one of probs and ageProbs}
}
\details{
We assume there are five stages for each patient: get infection (marked as Inf); get onset (marked as Ons); get medication (marked as Med); 
get recover (marked as Rec); go to hospitalization (marked as Hosp);  or eventually dead (marked as Death). There are also three levels of disease cases: mild, serious or deadly marked as M, S and D seperately.
}
\value{The value returned by \code{pandemicParams} is a list. Each list is for a different stage movement. Each list contains 4 values: mean, shape, zeros, scale. These values are chosen defaultly.

}
\references{
}
\author{
Patrick Brown
}
\note{

}


\seealso{ \link[addAgeProbs]{addAgeProbs}
}
\examples{
params = pandemicParams( 
  MedHospD = c(mean = 1.5, shape = 2, zeros = 0.2)
  )
params = pandemicParams(
  probs=NULL,
  ageProbs=list(
    S=addAgeProbs(),
    D=addAgeProbs()
  )
)
}

\keyword{ }
\keyword{ }% __ONLY ONE__ keyword per line
