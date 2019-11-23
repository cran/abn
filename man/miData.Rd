% mi-data.Rd ---
% Author           : Gilles Kratzer
% Created on :       17.07.2018
% Last modification :
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

\name{miData}
\alias{miData}
\title{Empirical Estimation of the Entropy from a Table of Counts}

\usage{
miData(freqs.table, method = c("mi.raw", "mi.raw.pc"))
}
\arguments{
  \item{freqs.table}{a table of counts.}
  \item{method}{a character determining if the Mutual Information should be normalized.}
  }
  
\description{
This function empirically estimates the Mutual Information from a table of counts using the observed frequencies.}

\details{
The mutual information estimation is computed from the observed frequencies through a plugin estimator based on entropy.

The plugin estimator is I(X, Y) = H (X) + H(Y) - H(X, Y), where H() is the entropy computed with \code{\link{entropyData}}.
}


\value{Mutual information estimate.}

\author{Gilles Kratzer}

\seealso{\code{\link{discretization}}}

\references{
Cover, Thomas M, and Joy A Thomas. (2012). "Elements of Information Theory". John Wiley & Sons.
}

\examples{
## Generate random variable
Y <- rnorm(n = 100, mean = 0, sd = 2)
X <- rnorm(n = 100, mean = 5, sd = 2)

dist <- list(Y="gaussian", X="gaussian")

miData(discretization(data.df = cbind(X,Y), data.dists = dist,
                discretization.method = "fd", nb.states = FALSE),
       method = "mi.raw")
}
\keyword{utilities}
