%                               -*- Mode: Rd -*- 
% simulateDag.Rd --- 
% Author          : Gilles Kratzer
% Created On      : 18.07.2018
% Last Modified By: 
% Last Modified On: 
% Update Count    : 
% Status          : Unknown, Use with caution!
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

\name{createDag}
\alias{createDag}
\alias{createAbnDag}
\encoding{latin1}
%- Also NEED an `\alias' for EACH other topic documented here.

\title{Create a legitimate DAGs}

\description{Create a legitimate DAG in the abn format.}

\usage{
createAbnDag( dag, data.df = NULL, data.dists = NULL, ...)
}
\arguments{
\item{dag}{a matrix or a formula specifying a DAG, see \sQuote{Details}.}
\item{data.df}{named dataframe or named vector.}
\item{data.dists}{named list giving the distribution for each node in the network. If not provided it will be sample and returned.}
\item{...}{further arguments passed to or from other methods.}
}

\details{An object of class \code{class(abnDag)} contains a named matrix describing the DAG and possibly additional objects such as the associated distributions of the nodes. 

If the dag is specified with a formula, either \code{data.df} or
\code{data.dists} is required with the \code{.} quantifier.

If the dag is specified with an unnamed matrix and both \code{data.df} and
\code{data.dists} are missing,  lower-case letters of the Roman alphabet
are used for the node names.
}

\value{An object of class \code{abnDag} containing a named matrix and a named list giving the distribution for each node.}


\examples{
createAbnDag( ~a+b|a, data.df=c("a"=1, "b"=1))

plot( createAbnDag( matrix( c(0,1,0,0),2,2))) 
}

\keyword{utilities}
\concept{DAG}
