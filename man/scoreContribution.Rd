%                               -*- Mode: Rd -*- 
% scoreContribution.Rd --- 
% Author          : Gilles Kratzer
% Created On      : 22/05/2019
% Last Modified By: 
% Last Modified On: 
% Update Count    : 
% Status          : Unknown, Use with caution!
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

\name{scoreContribution}
\alias{scoreContribution}
\encoding{latin1}

\title{Compute the score contribution of each individual observation in a network.}

\description{This function computes the score contribution of each individual observations to the total network score.}


\usage{
scoreContribution(object = NULL, 
                  dag.m = NULL, data.df = NULL, data.dists = NULL,
                  verbose = FALSE)
}
    

\arguments{
 \item{object}{an object of class '\code{abnLearned}' produced by \code{\link{mostprobable}}, \code{\link{search.heuristic}} or \code{\link{search.hillclimber}}.}
 \item{dag.m}{a matrix or a formula statement (see details) defining the network structure, a directed acyclic graph (DAG), see details for format. Note that colnames and rownames must be set.}
 \item{data.df}{a data frame containing the data used for learning the network, binary variables must be declared as factors and no missing values all allowed in any variable.}
 \item{data.dists}{a named list giving the distribution for each node in the network, see details.}
 \item{verbose}{if \code{TRUE} then provides some additional output.}
 
}

\details{
This function computes the score contribution of each individual observations to the total network score. This function is available only in the `mle` settings. To do so one uses the \code{\link{glm}} and \code{\link{predict}} functions. This funtion is an attempt to perform diagnostic for an ABN analysis.
}

\value{A named list that contains the scores controbutions: maximum likelihood, aic, bic, mdl and diagonal values of the hat matrix.}

\author{Gilles Kratzer}

\examples{

## Use a subset of a built-in simulated data set
mydat <- ex1.dag.data[,c("b1","g1","p1")]

## setup distribution list for each node
mydists <- list(b1="binomial", g1="gaussian", p1="poisson")

## now build cache
mycache <- buildscorecache(data.df = mydat, data.dists = mydists, max.parents = 2, method = "mle")

## Find the globally best DAG
mp.dag <- mostprobable(score.cache=mycache, score="bic", verbose = FALSE)

out <- scoreContribution(object = mp.dag)

## Observations contribution per network node
boxplot(out$bic)

}

\keyword{utilities}
\concept{abn}
