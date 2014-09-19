\name{corFDG}
\alias{corFDG}
\title{
   Sample pairwise dependence coefficients
}
\description{
  Computes sample (i.e., empirical) pairwise dependence
  coefficients such as the Kendall's tau, the Spearman's rho, and the
  upper tail dependence coefficient (the later being proved to have good
  statistical properties only for extreme-value copulas and in the case
  of known margins).
}
\usage{
corFDG(x, depcoefType)
}
\arguments{
  \item{x}{the data}
  \item{depcoefType}{type of dependence coefficient: one of
    'kendall', 'spearman', or 'utdc'}
}
\details{
  The function 'corFDG' simply calls 'cor()' when 'depcoefType' is
  either 'kendall' or 'spearman'. If it is 'utdc', then the sample upper
  tail dependence coefficient for extreme-value copulas introduced in
  the reference below is computed.
}
\value{
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
A matrix of size 'd' times 'd' where 'd' is the number of variables.
The element in the i-th row and j-th column is the dependence
coefficient between the i-th and j-th variable.
}
\references{
  Ferreira, M., Nonparametric estimation of the tail-dependence
  coefficient, REVSTAT--Statistical Journal, 2013
}
\author{
Gildas Mazo
}
%%\note{
%%  ~~further notes~~
%%}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
\code{\link{cor}}
}
\examples{
## Generate an extreme-value copula
myFDGcopula <- FDGcopula("frechet", c(.3,.5,.7,.9), extremevalue=TRUE)
dat <- rFDG(100, myFDGcopula)
## Compute the sample upper tail dependence coefficient given in the reference
corFDG(dat, "utdc")
}



