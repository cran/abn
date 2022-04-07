% build_control.Rd ---
% Author           : Kalina Cherneva
% Created on :       26/10/2021
% Last modification : 26/10/2021
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

\name{fit.control}

\alias{fit.control}

\title{
Control the iterations in \code{\link{fitAbn}}
}

\description{
  Allow the user to set restrictions in the \code{\link{fitAbn}} for both the
  Bayesian and the MLE approach.
}

\usage{
fit.control(method = "bayes", mean = 0, prec = 0.001, loggam.shape = 1,
            loggam.inv.scale = 5e-05, max.mode.error = 10, max.iters = 100,
            epsabs = 1e-07, error.verbose = FALSE, trace = 0L, epsabs.inner = 1e-06,
            max.iters.inner = 100, finite.step.size = 1e-07,
            hessian.params = c(1e-04, 0.01), max.iters.hessian = 10,
            max.hessian.error = 1e-04, factor.brent = 100, maxiters.hessian.brent = 10,
            num.intervals.brent = 100, min.pdf = 0.001, n.grid = 250, std.area = TRUE,
            marginal.quantiles = c(0.025, 0.25, 0.5, 0.75, 0.975), max.grid.iter = 1000,
            marginal.node = NULL, marginal.param = NULL, variate.vec = NULL,
            max.irls = 100, tol = 10^-11, seed = 9062019)
}

\arguments{
  \item{method}{a character that takes one of two values: "bayes" or "mle"}
  \item{mean}{the prior mean for all the Gaussian additive terms for each node.}
  \item{prec}{the prior precision for all the Gaussian additive terms for each node.}
  \item{loggam.shape}{the shape parameter in the Gamma distributed prior for the precision in any Gaussian nodes, also used for group-level precision is applicable.}
  \item{loggam.inv.scale}{the inverse scale parameter in the Gamma distributed prior for the precision in any Gaussian nodes, also used for group-level precision, is applicable. }
  \item{max.mode.error}{if the estimated modes from INLA differ by a factor of max.mode.error or more from those computed internally, then results from INLA are replaced by those computed internally. To force INLA always to be used, then \code{max.mode.error=100}, to force INLA never to be used \code{max.mod.error=0}. See details.}
  \item{max.iters}{total number of iterations allowed when estimating the modes in Laplace approximation}
  \item{epsabs}{absolute error when estimating the modes in Laplace approximation for models with no random effects.}
  \item{error.verbose}{logical, additional output in the case of errors occurring in the optimization}
  \item{trace}{Non-negative integer. If positive, tracing information on the progress of the "L-BFGS-B" optimization is produced. Higher values
          may produce more tracing information. (There are six levels of tracing.  To understand exactly what these do see the source code.)}
  \item{epsabs.inner}{absolute error in the maximization step in the (nested) Laplace approximation for each random effect term}
  \item{max.iters.inner}{total number of iterations in the maximization step in the nested Laplace approximation}
  \item{finite.step.size}{suggested step length used in finite difference estimation of the derivatives for the (outer) Laplace approximation when estimating modes}
  \item{hessian.params}{a numeric vector giving parameters for the adaptive algorithm, which determines the optimal step size in the finite-difference estimation of the Hessian. First entry is the initial guess, second entry absolute error}
  \item{max.iters.hessian}{integer, maximum number of iterations to use when determining an optimal finite difference approximation (Nelder-Mead)}
  \item{max.hessian.error}{if the estimated log marginal likelihood when using an adaptive 5pt finite-difference rule for the Hessian differs by more than max.hessian.error from when using an adaptive 3pt rule then continue to minimize the local error by switching to the Brent-Dekker root bracketing method, see details}
  \item{factor.brent}{if using Brent-Dekker root bracketing method then define the outer most interval end points as the best estimate of h (stepsize) from the Nelder-Mead as (h/factor.brent,h*factor.brent)}
  \item{maxiters.hessian.brent}{maximum number of iterations allowed in the Brent-Dekker method}
  \item{num.intervals.brent}{the number of initial different bracket segments to try in the Brent-Dekker method}
  \item{min.pdf}{the value of the posterior density function below which we stop the estimation only used when computing marginals, see details.}
  \item{n.grid}{recompute density on an equally spaced grid with \code{n.grid} points.}
  \item{std.area}{logical, should the area under the estimated posterior density be standardized to exactly one, useful for error checking.}
  \item{marginal.quantiles}{vector giving quantiles at which to compute the posterior marginal distribution at.}
  \item{max.grid.iter}{gives number of grid points to estimate posterior density at when not explicitly specifying a grid used to avoid excessively long computation.}
  \item{marginal.node}{used in conjunction with \code{marginal.param} to allow bespoke estimate of a marginal density over a specific grid. value from 1 to the number of nodes.}
 \item{marginal.param}{used in conjunction with \code{marginal.node}. value of 1 is for intercept, see modes entry in results for the appropriate number.}
  \item{variate.vec}{a vector containing the places to evaluate the posterior marginal density, must be supplied if \code{marginal.node} is not null}
  \item{max.irls}{integer given the maximum number of run for estimating network scores using an Iterative Reweighed Least Square algorithm.}
  \item{tol}{real number giving the minimal tolerance expected to terminate the Iterative Reweighed Least Square algorithm to estimate network score.}
  \item{seed}{a non-negative integer which sets the seed.}
}

\value{
A list with 26 components for the Bayesian approach, or a list with 3 components for "mle".
}


\examples{
ctrlmle <- fit.control(method = "mle", max.irls = 100, tol = 10^-11, seed = 9062019)

ctrlbayes <- fit.control(method = "bayes", mean = 0, prec = 0.001, loggam.shape = 1,
  loggam.inv.scale = 5e-05, max.mode.error = 10, max.iters = 100,
  epsabs = 1e-07, error.verbose = FALSE, epsabs.inner = 1e-06,
  max.iters.inner = 100, finite.step.size = 1e-07, hessian.params = c(1e-04, 0.01),
  max.iters.hessian = 10, max.hessian.error = 1e-04, factor.brent = 100,
  maxiters.hessian.brent = 10, num.intervals.brent = 100, min.pdf = 0.001,
  n.grid = 100, std.area = TRUE, marginal.quantiles = c(0.025, 0.25, 0.5, 0.75, 0.975),
  max.grid.iter = 1000, marginal.node = NULL, marginal.param = NULL, variate.vec = NULL,
  seed = 9062019)
}




