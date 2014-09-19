\name{ltdcFDG}
\alias{ltdcFDG}
\title{
  Lower tail dependence coefficient of FDG copulas 
}
\description{
  Calculates the lower tail dependence coefficient of FDG copulas.
}
\usage{
ltdcFDG(FDGcopula)
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
\code{\link{rhoFDG}}, \code{\link{tauFDG}}, \code{\link{utdcFDG}}
}
\examples{
## FDG copula with 'exponential' generators
myExpo <- FDGcopula("exponential", 1:4)
ltdcFDG(myExpo) 

## FDG copula with Cuadras-Auge generators
myCA <- FDGcopula("cuadrasauge", c(.3,.5,.7,.9))
ltdcFDG(myCA) # lower tail independent
}



