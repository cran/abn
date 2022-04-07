%                               -*- Mode: Rd -*-
% simulateAbn.Rd ---
% Author          : Gilles Kratzer
% Created On      : 17.01.2017
% Last Modified By:
% Last Modified On:
% Update Count    :
% Status          : Unknown, Use with caution!
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

\name{simulateAbn}
\alias{simulateAbn}
\alias{simulate.abn}
\alias{simulateabn}
\encoding{latin1}
%- Also NEED an `\alias' for EACH other topic documented here.

\title{Simulate from an ABN Network}

\description{Simulate one or more responses from an ABN network corresponding to a fitted object using a formula statement or an adjacency matrix.}

\usage{
simulateAbn(data.dists = NULL,
                       data.param = NULL,
                       data.param.var = NULL,
                       data.param.mult = NULL,
                       n.chains =10,
                       n.adapt = 1000,
                       n.thin = 100,
                       n.iter = 10000,
                       bug.file=NULL,
                       verbose=TRUE,
                       simulate=TRUE,
                       keep.file=FALSE,
                       seed=42)
}

%- maybe also `usage' for other objects documented here.
\arguments{
  \item{data.dists}{named list giving the distribution for each node in the network, see \sQuote{Details}.}
  \item{data.param}{named matrix, which have to be square with as many entries as the number of variables, each element is the coefficient (for specifications see \sQuote{Details}) used in the glm to simulate responses.}
  \item{data.param.var}{optional matrix, which should be square and having as many entries as number of variables, which contains the precision values for gaussian nodes. Default is set to 1.}
  \item{data.param.mult}{optional matrix, which should be square and having as many entries as number of variables, which contains the multinomial coefficient.}
  \item{n.chains}{number of parallel chains for the model.}
  \item{n.thin}{number of parallel chains for the model.}
  \item{n.iter}{number of iteration to monitor.}
  \item{n.adapt}{number of iteration for adaptation. If \code{n.adapt} is set to zero, then no adaptation takes place.}
  \item{bug.file}{name of the user specific bug file. If missing \code{"model.bug"} will be used. See \sQuote{Details}.}
  \item{verbose}{logical. Default \code{TRUE}. Should R report extra information on progress?}
  \item{simulate}{logical. Default \code{TRUE}. If set to \code{FALSE}, no simulation will be run only creation of the bug file.}
  \item{keep.file}{logical. Default \code{FALSE}. If set to \code{TRUE}, the bug file generated will be kept afterwards.}
  \item{seed}{by default set to 42.}
}

\details{
This function use \code{rjags} to simulate data from a DAG. It first creates a \code{bug} file, in the actual repository, then use it to simulate the data. This function output a data frame. The \code{bug} file can be run using \code{rjags} separately.

The coefficients given in the \code{data.param} are: the logit of the probabilities for binomial nodes, the means of the gaussian nodes, and the log of the Poison parameter. Additionally, a matrix \code{data.param.var} could give precision values for gaussian nodes (default is set to 1).

Binary and multinomial variables must be declared as factors, and the argument \code{data.dists} must be a list with named arguments, one for each of the variables in \code{data.df} (except a grouping variable - if present), where each entry is either \code{"poisson"}, \code{"binomial"}, \code{"multinomial"} or \code{"gaussian"}, see examples below. The Poisson distributed variables use log and Binomial and multinomial distributed ones the logit link functions. Note that "binomial" here actually means only binary, one Bernoulli trial per row in \code{data.df}.

The number of simulated data (rows of the outputted data frame) is given by \code{n.iter} divided by \code{n.thin}.

The bug file contains a description of the model in the JAGS dialect of the BUGS language.
It is possible to only generate this file and to reuse this for later or other purposes.
If a bug file name is specified and the file exists, it will be used.
If the file does not exist or \code{bug.file} is missing, the bug file
is created. Default name is \code{"model.bug"}.

The \code{verbose} argument is passed appropriately to the JAGS functions.

}

\value{A data frame containing simulated data.}

\references{
Further information about \pkg{abn} can be found at:\cr
  \url{http://r-bayesian-networks.org}}

\author{Gilles Kratzer}


\examples{
## Define set of distributions:
dist<-list(a="gaussian", b="gaussian", c="gaussian", d="gaussian",
           e="binomial", f="binomial")

## Define parameter matrix:
data.param <- matrix(c(1,2,0.5,0,20,0,
                       0,1,3,10,0, 0.8,
                       0,0,1,0,0,0,
                       0,0,0,1,0,0,
                       0,0,0,0,0.5,1,
                       0,0,0,0,0,0), nrow = 6L, ncol = 6L, byrow = TRUE)

## Define precision matrix:
data.param.var <- matrix(0, nrow = 6L, ncol = 6L)
diag(data.param.var) <- c(10,20,30,40,0,0)

## Plot the dag
plotAbn(dag = ~a|b:c:e+b|c:d:f+e|f, data.dists = dist)

## Simulate the data
out <- simulateAbn(data.dists=dist, n.chains=1, n.thin=1, n.iter=1000,
            data.param=data.param, data.param.var=data.param.var)
}

\keyword{utilities}
\concept{DAG}
