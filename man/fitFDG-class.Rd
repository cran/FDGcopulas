\name{fitFDG-class}
\Rdversion{1.1}
\docType{class}
\alias{fitFDG-class}

\title{Class \code{"fitFDG"}}
\description{
  Class to represent fitted FDG copulas
  %%  ~~ A concise (1-5 lines) description of what the class is. ~~
}
\section{Slots}{
  \describe{
    \item{\code{estimate}:}{Object of class \code{"numeric"} estimated
      parameter vector }
    \item{\code{var.est}:}{Object of class \code{"matrix"} estimated
      variance-covariance matrix of scaled estimator }
    \item{\code{optimalvalues}:}{Object of class \code{"numeric"}
      optimal values of the loss function which was minimized during the
      estimation procedure}
    \item{\code{convergence}:}{Object of class \code{"list"} monitoring
      parameters returned by 'optim()'}
    \item{\code{FDGcopula}:}{Object of class \code{"FDGcopula"}
      estimated copula}
  }
}
\section{Methods}{
  No methods defined with class "fitFDG" in the signature.
}
\section{Details}{
  The estimated variance-covariance matrix is that of the estimator
  times the square root of the sample size.
}
\references{
  [1] Mazo G., Girard, S., Forbes, F. A flexible and tractable class of
  one-factor copulas, http://hal.archives-ouvertes.fr/hal-00979147
  
  [2] Mazo G., Girard, S., Forbes, F. Weighted least-squares inference based
  on dependence coefficients for multivariate copulas,
  http://hal.archives-ouvertes.fr/hal-00979151
}
\author{
  Gildas Mazo
}
%\note{
%%  ~~further notes~~
%}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
  \code{\link{fitFDG}} for fitting FDG copulas
}
\examples{
showClass("fitFDG")
}
\keyword{classes}
