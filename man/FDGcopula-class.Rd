\name{FDGcopula-class}
\Rdversion{1.1}
\docType{class}
\alias{FDGcopula-class}

\title{Class \code{"FDGcopula"}}
\description{
Class to represent a FDG copula
}
\section{Slots}{
  \describe{
    \item{\code{dimension}:}{Object of class \code{"integer"} Dimension
      of the copula}
    \item{\code{parameters}:}{Object of class \code{"numeric"} Vector of
    parameters}
    \item{\code{family}:}{Object of class \code{"character"} Family of
      the generators}
    \item{\code{extremevalue}:}{Object of class \code{"logical"} Boolean
      indicating if the copula is an extreme-value copula}
    \item{\code{parameterrange}:}{Object of class \code{"numeric"}
      Interval containing the bounds of the parameters}
  }
}
\section{Methods}{
No methods defined with class "FDGcopula" in the signature.
}
\references{
  Mazo G., Girard, S., Forbes, F. A flexible and tractable class of
  one-factor copulas, http://hal.archives-ouvertes.fr/hal-00979147
}
\author{
Gildas Mazo
}

\seealso{
\code{\link{FDGcopula}} for constructing FDG copulas
}
\examples{
showClass("FDGcopula")
}
\keyword{classes}
