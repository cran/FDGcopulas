\name{fitFDG}
\alias{fitFDG}
\title{
  Estimation of FDG copulas
}
\description{
  Estimates FDG copulas.
}
\usage{
fitFDG(FDGcopula, data, depcoefType="spearman", nbInit=1,
                   W=NA, method="L-BFGS-B", estimate.variance=TRUE,
                   nb.rep=100, nb.obs=100, dcData=NA, sizeSubSample=10000)
}
\arguments{
  \item{FDGcopula}{copula to fit the data}
  \item{data}{data to be fitted}
  \item{depcoefType}{type of dependence coefficient to be used in the
    estimation method}
  \item{nbInit}{number of initialisations in the optimization algorithm
    to find the global optimum}
  \item{W}{weight matrix to be used in the estimation method}
  \item{method}{the 'method' for 'optim()'}
  \item{estimate.variance}{boolean indicating if the asymptotic
    variance-covariance matrix should be estimated}
  \item{nb.rep}{number of replications to be used in the estimation of the
    asymptotic variance-covariance matrix (it has no effect if
    estimate.variance=FALSE)}
  \item{nb.obs}{size of the simulated samples on which is based the
    estimation of the asymptotic variance-covariance matrix}
  \item{dcData}{if not NA, matrix of size 'd' times 'd', where 'd' is
    the dimension of the fitted copula, consisting of the sample
    pairwise dependence coefficients that the user wants to use for the
    estimation method}
  \item{sizeSubSample}{size of the sample over which is to be
    taken the maximum when generating extreme-value copulas (needed if
    estimate.variance=TRUE and FDGcopula@extremevalue=TRUE,
    no effect otherwise)}
}
\details{
  The method used to estimate the parameters of FDG copulas is a
  weighted least squares estimator based on dependence coefficients (see
  [2]). The
  coefficients implemented are the Spearman's rho, the Kendall's tau,
  and the upper tail dependence coefficient in the case of extreme-value
  copulas. If the user wishes to use other coefficients, it is possible but he/she
  should provide his own sample pairwise dependence coefficients with the
  matrix \code{dcData}. The estimation of the asymptotic
  variance-covariance matrix of 'sqrt(n)(theta hat - theta)', where 'n' is
  the sample size, 'theta' is the parameter vector, and 'theta hat' is the
  weighted least square estimator, is carried out by simulation. More
  precisely, \code{nb.rep} replications of datasets of size \code{nb.obs} are
  simulated according to the fitted FDG copula. For each dataset,
  the sample dependence coefficients are calculated, and, then, their sample variances / covariances are computed. In the case where the upper tail
  dependence coefficients were chosen to perform the estimation, a
  different approximation is used. Since the
  margins are assumed to be known, there is a simple formula
  for the variances / covariances given in (15) of [2]. These quantities
  within this formula can be approximated by standard empirical means
  calculated on a single big dataset from the underlying extreme-value
  copula. To simulate that dataset, the variable \code{sizeSubSample} is used
  along with \code{nb.rep}: \code{nb.rep} sub-samples of size \code{sizeSubSample} are
  simulated, and for each sub-sample, the maximum is taken, thus leading
  to a final dataset of size \code{nb.rep}. The empirical means to
  approximate the asymptotic variances / covariances are computed on
  this last final dataset.
}
\value{
A \code{fitFDG} class object containing the slots:
\item{estimate}{the estimated parameter vector}
\item{var.est}{the asymptotic variance-covariance matrix}
\item{optimalValues}{the optimal value(s) of the loss function}
\item{convergence}{the output monitor parameters returned by 'optim()'}
\item{copula}{an object of the same class as FDGcopula, where the
  slot containing the parameters is filled with the estimated parameters}
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
\seealso{\code{\link{fitFDG-class}}
}
\examples{
## Create an object of class 'FDGcopula'
theta <- c(.3,.5,.7,.9)
myFDGcopula <- FDGcopula("frechet", theta)

## Generate a sample from a FDG copula with Frechet generators
## and parameter vector 'theta'
dat <- rFDG(100, myFDGcopula)
## Fit a FDG copula to the data
myFittedCopula <- fitFDG(myFDGcopula, dat)
myFittedCopula
}



