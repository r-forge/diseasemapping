\name{simLGCP}
\alias{simLGCP}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Simulate a log-Gaussian Cox process
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
simLGCP(N = 101, mu = 3, beta = c(1, 1.5), sigma = 1, range = 0.3, roughness = 2)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{N}{
Number of grid cells in the X and Y directions on which the intensity surface is simulated
}
  \item{mu}{
Intercept parameter (on the log scale)
}
  \item{beta}{
Coefficients for two spatial covariates (see details)
}
  \item{sigma}{
Standard deviation of the latent spatial process
}
  \item{range}{
Range parameter for the latent spatial process
}
  \item{roughness}{
Matern roughness parameter for the latent spatial process
}

}
\details{
The intensity surface is $O(s)\lambda(s)$ with
\[
\lambda(s) = \exp[\mu + \beta_1 X_1(s) + \beta_2 X_2(s) + U(s)]
\]
and $U$ being a Gaussian random field with
\[
cov[U(s+h), U(s)] = \sigma^2 M(||h||/\phi;\nu)
\]
where $\nu$ is the roughness parameter and $\phi$ is the range parameter.
}
\value{
 \item{points }{A \code{SpatialPoints} object containing event locations.}
 \item{params }{A list of model parameters}
 And the raster objects
 \item{lambda }{The relative intensity $\lambda(s)$}
 \item{intensity}{The intensity $O(s)\lambda(s)$}
 \item{U}{The latent Gaussian surface}
 \item{covs}{The two spatial covariates}
 \item{offset}{The offset surface $O(s)$}
}
\examples{

mysim = simLGCP(N=20)
image(mysim$intensity)
points(mysim$points)
 
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{simulation}
