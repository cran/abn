% build_control.Rd ---
% Author           : Kalina Cherneva
% Created on :       04/10/2021
% Last modification : 04/10/2021
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

\name{build.control}

\alias{build.control}

\title{
Control the iterations in \code{\link{buildScoreCache}}
}

\description{
  Allow the user to set restrictions in the \code{\link{buildscorecache}} for both the
  Bayesian and the MLE approach.
}

\usage{
build.control(method = "bayes", max.mode.error = 10, mean = 0, prec = 0.001,
          loggam.shape = 1, loggam.inv.scale = 5e-05, max.iters = 100, epsabs = 1e-07,
          error.verbose = FALSE, trace = 0L, epsabs.inner = 1e-06, max.iters.inner = 100,
          finite.step.size = 1e-07, hessian.params = c(1e-04, 0.01),
          max.iters.hessian = 10, max.hessian.error = 0.5, factor.brent = 100,
          maxiters.hessian.brent = 100, num.intervals.brent = 100,
          ncores = 0, max.irls = 100, tol = 10^-8, seed = 9062019)
}

\arguments{
  \item{method}{a character that takes one of two values: "bayes" or "mle"}
  \item{max.mode.error}{if the estimated modes from INLA differ by a factor of max.mode.error or more from those computed internally, then results from INLA are replaced by those computed internally. To force INLA always to be used, then max.mode.error=100, to force INLA never to be used max.mod.error=0.}
  \item{mean}{the prior mean for all the Gaussian additive terms for each node}
  \item{prec}{the prior precision for all the Gaussian additive term for each node}
  \item{loggam.shape}{the shape parameter in the Gamma distribution prior for the precision in a Gaussian node}
  \item{loggam.inv.scale}{the inverse scale parameter in the Gamma distribution prior for the precision in a Gaussian node}
  \item{max.iters}{total number of iterations allowed when estimating the modes in Laplace approximation}
  \item{epsabs}{absolute error when estimating the modes in Laplace approximation for models with no random effects.}
  \item{error.verbose}{logical, additional output in the case of errors occurring in the optimization}
  \item{trace}{Non-negative integer. If positive, tracing information on the progress of the "L-BFGS-B" optimization is produced. Higher values
          may produce more tracing information. (There are six levels of tracing.  To understand exactly what these do see the source code.)}
  \item{epsabs.inner}{absolute error in the maximization step in the (nested) Laplace approximation for each random effect term}
  \item{max.iters.inner}{total number of iterations in the maximization step in the nested Laplace approximation}
  \item{finite.step.size}{suggested step length used in finite difference estimation of the derivatives for the (outer) Laplace approximation when estimating modes}
  \item{hessian.params}{a numeric vector giving parameters for the adaptive algorithm, which determines the optimal stepsize in the finite-difference estimation of the hessian. First entry is the initial guess, second entry absolute error}
  \item{max.iters.hessian}{integer, maximum number of iterations to use when determining an optimal finite difference approximation (Nelder-Mead)}
  \item{max.hessian.error}{if the estimated log marginal likelihood when using an adaptive 5pt finite-difference rule for the Hessian differs by more than max.hessian.error from when using an adaptive 3pt rule then continue to minimize the local error by switching to the Brent-Dekker root bracketing method}
  \item{factor.brent}{if using Brent-Dekker root bracketing method then define the outer most interval end points as the best estimate of h (stepsize) from the Nelder-Mead as (h/factor.brent,h*factor.brent)}
  \item{maxiters.hessian.brent}{maximum number of iterations allowed in the Brent-Dekker method}
  \item{num.intervals.brent}{the number of initial different bracket segments to try in the Brent-Dekker method}
  \item{max.irls}{total number of iterations for estimating network scores using an Iterative Reweighed Least Square algorithm}
  \item{tol}{real number giving the minimal tolerance expected to terminate the Iterative Reweighed Least Square algorithm to estimate network score.}
  \item{ncores}{The number of cores to parallelize to, see \sQuote{Details}.}
  \item{seed}{a non-negative integer which sets the seed.}
}
\details{
  Parallelization over all children is possible via the function \code{foreach} of the package \pkg{doParallel}. \code{ncode=1} uses single threaded \code{foreach}. \code{ncode=-1} uses all available cores but one.

  With \code{ncores=0} a simple \code{for} loop is used.

}
\value{
A list with 18 components for the Bayesian approach, or a list with 4 components for "mle"
}

\examples{
ctrlmle <- build.control(method = "mle", ncores = 0, max.irls = 100, tol = 10^-11, seed = 9062019)

ctrlbayes <- build.control(method = "bayes", max.mode.error = 10, mean = 0, prec = 0.001,
   loggam.shape = 1, loggam.inv.scale = 5e-05, max.iters = 100,
   epsabs = 1e-07, error.verbose = FALSE, epsabs.inner = 1e-06,
   max.iters.inner = 100, finite.step.size = 1e-07, hessian.params = c(1e-04, 0.01),
   max.iters.hessian = 10, max.hessian.error = 0.5, factor.brent = 100,
   maxiters.hessian.brent = 100, num.intervals.brent = 100,
   tol = 10^-8, seed = 9062019)


}




