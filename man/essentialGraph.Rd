%                               -*- Mode: Rd -*-
% essentialGraph.Rd ---
% Author          : Gilles Kratzer
% Created On      : 17/07/2018
% Last Modified By: GK (rd file + naming + examples)
% Last Modified On:
% Update Count    :
% Status          : Unknown, Use with caution!
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

\name{essentialGraph}
\alias{essentialGraph}
\encoding{latin1}
%- Also NEED an `\alias' for EACH other topic documented here.

\title{Construct the essential graph}

\description{Constructs different versions of the essential graph from a given DAG}

\usage{essentialGraph(dag, node.names = NULL, PDAG = "minimal")
}

%- maybe also `usage' for other objects documented here.
\arguments{
  \item{dag}{a matrix or a formula statement (see \sQuote{Details} for format) defining the network structure, a directed acyclic graph (DAG).}
  \item{node.names}{a vector of names if the DAG is given via formula, see \sQuote{Details}.}
  \item{PDAG}{a character value that can be: minimal or complete, see \sQuote{Details}.}
}


\details{

This function returns an essential graph from a DAG, aka acyclic partially directed graph (PDAG). This can be useful if the learning procedure is defined up to a Markov class of equivalence. A minimal PDAG is defined as only directed edges are those who participate in v-structure. Whereas the completed PDAG: every directed edge corresponds to a compelled edge, and every undirected edge corresponds to a reversible edge.

The \code{dag} can be provided using a formula statement (similar to glm). A typical formula is \code{ ~ node1|parent1:parent2 + node2:node3|parent3}. The formula statement have to start with \code{~}. In this example, node1 has two parents (parent1 and parent2). node2 and node3 have the same parent3. The parents names have to exactly match those given in \code{node.names}. \code{:} is the separator between either children or parents, \code{|} separates children (left side) and parents (right side), \code{+} separates terms, \code{.} replaces all the variables in \code{node.names}.
}

\value{A matrix giving the PDAG.}

\references{
West, D. B. (2001). Introduction to Graph Theory. Vol. 2. Upper Saddle River: Prentice Hall.

Chickering, D. M. (2013) A Transformational Characterization of Equivalent Bayesian Network Structures, arXiv:1302.4938.



Further information about \pkg{abn} can be found at:\cr
  \url{http://r-bayesian-networks.org}}

\author{Gilles Kratzer}


\examples{
dag <- matrix(c(0,0,0, 1,0,0, 1,1,0), nrow = 3, ncol = 3)
dist <- list(a="gaussian", b="gaussian", c="gaussian")
colnames(dag) <- rownames(dag) <- names(dist)

essentialGraph(dag)
}
\keyword{utilities}
\concept{DAG}
