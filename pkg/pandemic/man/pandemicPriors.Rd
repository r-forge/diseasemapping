\name{pandemicPriors}
\Rdversion{1.1}
\alias{pandemicPriors}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Get priors  
}
\description{
Getting the priors of different disease stages for the MCMC simulation. 
}
\usage{
pandemicPriors(InfOns = meanShapeZerosPrior(), OnsMedM = meanShapeZerosPrior(), OnsMedS = meanShapeZerosPrior(), OnsMedD = meanShapeZerosPrior(), MedRec = meanShapeZerosLostPrior(), MedHospS = meanShapeZerosPrior(), MedHospD = meanShapeZerosPrior(), HospRec = meanShapeZerosPrior(), HospDeath = meanShapeZerosPrior(), probs = probsPrior())
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{InfOns}{prior for time from stage infection to onsets. }
  \item{OnsMedM}{prior for time from stage onsets to on medication for mild cases.}
  \item{OnsMedS}{prior for time from onsets to on medication for serious cases.}
  \item{OnsMedD}{prior for time from onsets to on medication for deadly cases.}
  \item{MedRec}{prior for time from on medication to recover.}
  \item{MedHospS}{prior for time from on medication to hospitalization for serious cases.}
  \item{MedHospD}{prior for time from on medication to hospitalization for deadly cases.}
  \item{HospRec}{prior for time from hospitalization to recover.}
  \item{HospDeath}{prior for time from hospitalization to death.}
  \item{probs}{prior for probabily of mild cases, serious cases and deadly cases.}
}
\details{
We assume there are five stages for each patient: get infection (marked as Inf); get onset (marked as Ons); get medication (marked as Med); 
get recover (marked as Rec); go to hospitalization (marked as Hosp);  or eventually dead (marked as Death). There are also three levels of disease cases: mild, serious or deadly marked as M, S and D seperately. This is a function used to get the prior distributions for the time between different disease stages.
}
\value{
The function returns a list. 
}
\references{

}
\author{
Patrick Brown
}
\note{

}



\seealso{

}
\examples{

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{}
\keyword{}% __ONLY ONE__ keyword per line
