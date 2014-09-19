\name{rFDG}
\alias{rFDG}
\title{
  Simulation of FDG copula models
}
\description{
  Simulates from FDG copula models.
}
\usage{
rFDG(n, FDGcopula, sizeSubSample=10000)
}
\arguments{
  \item{n}{size of the sample to be generated}
  \item{FDGcopula}{copula from which the sample is to be generated}
  \item{sizeSubSample}{
    size of the sub-samples over which is to be taken the maximum when
    generating extreme-value copulas
  }
}
\details{
  An observation of an extreme-value distribution is essentially the
  maximum over a sub-sample of (theoretically) infinite size. In practice,
  the size of this sub-sample is set to \code{sizeSubSample}. 
}
\value{
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
A matrix of size 'n' times the dimension of the copula.
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

%\seealso{
%\code{\link{FDGcopula-class}}
%}
\examples{
## Generate a sample from a FDG copula with Fr\'echet generators
## and parameter vector 'theta'
theta <- c(.3,.5,.7,.9)
myFDGcopula <- FDGcopula("frechet", theta)
dat <- rFDG(30, myFDGcopula)
head(dat)
}



