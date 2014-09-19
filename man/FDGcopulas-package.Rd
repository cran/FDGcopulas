\name{FDGcopulas-package}
\alias{FDGcopulas-package}
\docType{package}
\title{
Deals with FDG copulas
}
\description{
Constructs, simulates from, calculates dependence
coefficients, and fits FDG copulas.
}
\details{
\tabular{ll}{
Package: \tab FDGcopulas\cr
Type: \tab Package\cr
Version: \tab 1.0\cr
Date: \tab 2014-09-19\cr
License: \tab GPL >=3 \cr
}
See the examples below to have an overview of how what the package can offer.
}
\author{
  Gildas Mazo
  Stephane Girard
  Maintainer: Gildas Mazo <gildas.mazo@free.fr>
}
\note{Feel free to contact the authors if have any comments, suggestions, or want to report a bug.}
\references{
  Mazo G., Girard, S., Forbes, F. A flexible and tractable class of
  one-factor copulas, http://hal.archives-ouvertes.fr/hal-00979147
}
\keyword{copula}
\keyword{multivariate}
%% \seealso{
%% ~~ Optional links to other man pages, e.g. ~~
%% ~~ \code{\link[<pkg>:<pkg>-package]{<pkg>}} ~~
%% }
\examples{
## creates an object of class 'FDGcopula' 
myFDGcopula <- FDGcopula("frechet", c(.3,.5,.7,.9))

## compute the pairwise dependence coefficients 
## Spearman's rho:
rhoFDG(myFDGcopula)
## Kendall's tau:
tauFDG(myFDGcopula)
## Upper tail dependence coefficient:
utdcFDG(myFDGcopula)
## Lower tail dependence coefficient:
ltdcFDG(myFDGcopula)

## simulates data ##
dat <- rFDG(30, myFDGcopula)

## fit data ##
myFittedCopula <- fitFDG(myFDGcopula, dat)
}
