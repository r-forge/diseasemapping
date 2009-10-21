\name{ECS}
\Rdversion{1.1}
\alias{ECS}
\docType{data}
\title{
ECS data
}
\description{
Data from the Ontario portion of the ECS case-control study.  
}
\usage{data(ECS)}
\format{
  A data frame with 4396 observations on the following 17 variables.
  \describe{
      \item{\code{id}}{individual identifier.}
    \item{\code{date}}{a character vector of date the questionnaire was answered.}
    \item{\code{diagage}}{a character vector of diagnosis age. The default diagnosis date is July 1, 1996.}
    \item{\code{sex}}{a character vector of sex, "M" or "F".}
    \item{\code{smk_ever}}{TRUE if the individual is a current or former smoker, FALSE otherwise}
    \item{\code{smk_age}}{age the person started smoking.}
    \item{\code{smk_yrs}}{reported number of years smoked.}
    \item{\code{smk_now}}{TRUE if the individual is a current or former smoker, FALSE otherwise.}
    \item{\code{smk_avg}}{average number of cigarettes per day over entire time smoking.}
    \item{\code{smk_quit}}{the age the person quitted smoking.}
    \item{\code{packyrs}}{number of packs per year. Using 25 cigarettes as pack size. }
    \item{\code{post_code}}{a character vector of where they were living when they got the questionnaire.}
    \item{\code{EA1996}}{an integer vector of EA cecus tract from 1996.} 
    \item{\code{LAT}}{an integer vector of lattitude. }
    \item{\code{LONG}}{an integer vector of longitude.}
    \item{\code{CD}}{census division code.}
    \item{\code{cancer}}{a factor of the cancer type, with levels \code{control}  \code{stomach} \code{lung} \code{breast}}
    \item{\code{smoke}}{a factor of smoking status with levels \code{current} \code{former} \code{never}}
}
}
\details{
}
\source{
}
\references{
http://www.phac-aspc.gc.ca/publicat/cdic-mcc/21-1/g_e.html
}
\examples{
data(ECS)
head(ECS)

lungdata = ECS[ECS$cancer \%in\% c("lung","control"),]
lungdata$y = (lungdata$cancer == "lung")
summary(glm(y ~ smoke + diagage, data=lungdata))
}
\keyword{datasets}
