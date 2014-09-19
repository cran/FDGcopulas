\name{FDGcopula}
\alias{FDGcopula}
\title{
  Construction of FDG copula class object 
}
\description{
  Constructs a FDG copula class object.
}
\usage{
FDGcopula(family, parameters, extremevalue=FALSE, checkbounds=TRUE)
}
\arguments{
  \item{family}{character to indicate the family of the generators}
  \item{parameters}{the parameter vector}
  \item{extremevalue}{boolean indicating if the extreme-value limit
    copula corresponding to the FDG copula is
    to be constructed instead}
  \item{checkbounds}{boolean indicating if an error message is to be
    produced in case of parameter values are out of their theoretical bounds}
}
\details{
  The possible families are 'cuadrasauge', 'frechet', 'sinus' and 'exponential'.  Control over the possibility of NOT returning an error message in case
  of wrong parameters is given to be able to use certain optimization algorithms.}
\value{
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
An object of class \code{\link{FDGcopula-class}}
}
\references{
  Mazo G., Girard, S., Forbes, F. A flexible and tractable class of
  one-factor copulas, http://hal.archives-ouvertes.fr/hal-00979147
}
\author{
Gildas Mazo
}
%%\note{
%%  ~~further notes~~
%%}

%% ~Make other sections like Warning with \section{Warning }{....} ~
\seealso{
  \code{\link{FDGcopula-class}} for a description of the class
}
\examples{
## Create a FDG copula with Cuadras-Aug\'e generators
## and parameter vector 'theta'
theta <- c(.3,.5,.7,.9)
myFDGcopula <- FDGcopula("cuadrasauge", theta)
str(myFDGcopula)

## Create its extreme-value limit copula
myExtremeValueFDGcopula <- FDGcopula("cuadrasauge", theta,
extremevalue=TRUE)
str(myExtremeValueFDGcopula)

## Note: the two above copulas are NOT the same!
}



