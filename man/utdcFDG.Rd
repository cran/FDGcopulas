\name{utdcFDG}
\alias{utdcFDG}
\title{
  Upper tail dependence coefficient of FDG copulas 
}
\description{
  Calculates the upper tail dependence coefficient of FDG copulas.
}
\usage{
utdcFDG(FDGcopula)
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
\code{\link{rhoFDG}}, \code{\link{tauFDG}}, \code{\link{ltdcFDG}}
}
\examples{
## FDG copula with 'exponential' generators
myFDGcopula <- FDGcopula("exponential", 1:4)
utdcFDG(myFDGcopula) # upper tail independent
## FDG copula with 'sinus' generators
mySinus <- FDGcopula("sinus", c(.9,1.2,1.3,1.55))
utdcFDG(mySinus)
}


