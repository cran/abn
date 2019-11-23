%                               -*- Mode: Rd -*- 
% link_strength.Rd --- 
% Author          : Gilles Kratzer
% Created On      : 17/07/2018
% Last Modified By: 
% Last Modified On: 
% Update Count    : 
% Status          : Unknown, Use with caution!
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

\name{linkStrength}
\alias{linkStrength}
\alias{link.strength}
\alias{linkstrength}
\encoding{latin1}
%- Also NEED an `\alias' for EACH other topic documented here.

\title{A function that returns the strengths of the edge connections in a Bayesian Network learned from observational data.}

\description{A flexible implementation of multiple proxy for strength measures useful for visualizing the edge connections in a Bayesian Network learned from observational data.}

\usage{
linkStrength(dag, data.df = NULL, data.dists = NULL, 
              method = c("mi.raw", "mi.raw.pc", "mi.corr", "ls", "ls.pc", "stat.dist"),
              discretization.method = "doane")
}
%- maybe also `usage' for other objects documented here.
\arguments{
  \item{dag}{a matrix or a formula statement (see details for format) defining the network structure, a directed acyclic graph (DAG). Note that rownames must be set or given in \code{data.dist} if the DAG is given via a formula statement.}
  \item{data.df}{a data frame containing the data used for learning each node, binary variables must be declared as factors.}
  \item{data.dists}{a named list giving the distribution for each node in the network, see \sQuote{Details}.}
  \item{method}{the method to be used. See \sQuote{Details}.}
  \item{discretization.method}{a character vector giving the discretization method to use. See \code{\link{discretization}}.}
}


\details{This function returns multiple proxies for estimating the connection strength of the edges of a possibly discretized Bayesian network's dataset. The retuned connection strength measures are: the Raw Mutual Information (\code{mi.raw}), the Percentage Mutual information (\code{mi.raw.pc}), the Raw Mutual Information computed via correlation (\code{mi.corr}), the link strength (\code{ls}), the percentage link strength (\code{ls.pc}) and the statistical distance (\code{stat.dist}).

The general concept of entropy is defined for probability distributions. The probability is estimated from data using frequency tables. Then the estimates are plug-in in the definition of the entropy to return the so-called empirical entropy. A standard known problem of empirical entropy is that the estimations are biased due to the sampling noise. This is also known that the bias will decrease as the sample size increases.
The mutual information estimation is computed from the observed frequencies through a plug-in estimator based on entropy. For the case of an arc going from the node X to the node Y and the remaining set of parent of Y is denoted as Z.

The mutual information is defined as I(X, Y) = H(X) + H(Y) - H(X, Y), where H() is the entropy.

The Percentage Mutual information is defined as PI(X,Y) = I(X,Y)/H(Y|Z).

The Mutual Information computed via correlation is defined as MI(X,Y) = -0.5 log(1-cor(X,Y)^2).

The link strength is defined as LS(X->Y) = H(Y|Z)-H(Y|X,Z).

The percentage link strength is defined as PLS(X->Y) = LS(X->Y) / H(Y|Z).

The statistical distance is defined as SD(X,Y) = 1- MI(X,Y) / max(H(X),H(Y)).
}

\value{The function returns a named matrix with the requested metric.}

\references{
Boerlage, B. (1992).  Link strength in Bayesian networks. Diss. University of British Columbia.

Ebert-Uphoff, Imme. "Tutorial on how to measure link strengths in discrete Bayesian networks." (2009).

Further information about \pkg{abn} can be found at:\cr
  \url{http://r-bayesian-networks.org}}


\author{Gilles Kratzer}


\examples{

dist <- list(a="gaussian", b="gaussian", c="gaussian")
data.param <- matrix(c(0,1,0, 0,0,1, 0,0,0), nrow = 3L, ncol = 3L, byrow = TRUE)
    
data.param.var <- matrix(0, nrow = 3L, ncol = 3L)
diag(data.param.var) <- c(0.1,0.1,0.1)
    
out <- simulateAbn(data.dists = dist,
    n.chains = 1, n.adapt = 1000, n.thin = 1, n.iter = 100,
    data.param = data.param, data.param.var = data.param.var)

linkStrength(data.param, data.df = out, data.dists = dist,
             method = "ls", discretization.method = "sturges")

}

\keyword{utilities}
\concept{DAG}
