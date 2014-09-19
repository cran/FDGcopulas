\name{rhoFDG}
\alias{rhoFDG}
\title{
  Spearman's rho of FDG copulas 
}
\description{
  Calculates the Spearman's rho of FDG copulas.
}
\usage{
rhoFDG(FDGcopula)
}
\arguments{
  \item{FDGcopula}{the FDG copula class object}
}
\value{
A matrix of size 'd' times 'd', where 'd' is the dimension of the
copula. The element in the i-th row and j-th column is the dependence
coefficient of the i-th and j-th variable.
}
\references{
  Mazo G., Girard, S., Forbes, F. A flexible and tractable class of
  one-factor copulas, http://hal.archives-ouvertes.fr/hal-00979147
}
\author{
Gildas Mazo
}

\seealso{
\code{\link{tauFDG}}, \code{\link{utdcFDG}}, \code{\link{ltdcFDG}}
}
\examples{
## FDG copula with Frechet generators
myFDGcopula <- FDGcopula("frechet", c(.3,.5,.7,.9))
rhoFDG(myFDGcopula)
}



