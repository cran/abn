%                               -*- Mode: Rd -*- 
% simulateDag.Rd --- 
% Author          : Gilles Kratzer
% Created On      : 18.07.2018
% Last Modified By: 
% Last Modified On: 
% Update Count    : 
% Status          : Unknown, Use with caution!
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

\name{simulateDag}
\alias{simulateDag}
\encoding{latin1}
%- Also NEED an `\alias' for EACH other topic documented here.

\title{Simulate DAGs}

\description{Simulate a Directed Acyclic Graph (ABN) with arbitrary arc density.}

\usage{
simulateDag(node.name = NULL, 
            data.dists = NULL, 
            nn = 0.5)
}
    
%- maybe also `usage' for other objects documented here.
\arguments{
\item{node.name}{a vector of character giving the names of the nodes. It gives the size of the simulated DAG.}
\item{data.dists}{named list giving the distribution for each node in the network. If not provided it will be sample and returned.}
\item{nn}{a real number between 0 and 1 giving the network density.}
}

\details{
This function generates DAGs by sampling triangular matrices and reorder randomly columns and rows. The network density (\code{nn}) is used column wise as binomial sampling probability. Then the matrix is named using the user provided names.}

\value{A list containing a named matrix and a named list giving the distribution for each node.}

\references{
Further information about \bold{abn} can be found at:\cr
  \url{http://www.r-bayesian-networks.org}}

\author{Gilles Kratzer}

\examples{
##example using Ozon entries:

dist <- list(Ozone="gaussian",
Solar.R="gaussian",
Wind="gaussian",
Temp="gaussian",
Month="gaussian",
Day="gaussian")
simulateDag(node.name = names(dist),data.dists = dist,nn = 0.8)
}

\keyword{simulated structure, DAG}