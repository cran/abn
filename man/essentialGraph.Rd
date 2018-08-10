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

\title{Plot an ABN graphic}

\description{Plot an ABN DAG using formula statement or a matrix in using Rgraphviz through the graphAM class}

\usage{essentialGraph(dag.m=NULL, 
                      name=NULL, 
                      PDAG = "minimal")
}
    
%- maybe also `usage' for other objects documented here.
\arguments{
  \item{dag.m}{a matrix or a formula statement (see details for format) defining the network structure, a directed acyclic graph (DAG).}
  \item{name}{a vector of names if the DAG is given via formula, see details.}
  \item{PDAG}{a character value that can be: minimal or complete, see details.}
}


\details{

This function returns an essential graph from a DAG. This can be useful if the learning procedure is defined up to a Markov class of equivalence. A minimal PDAG is defined as only directed edges are those who participate in v-structure. Whereas the completed PDAG: every directed edge corresponds to a compelled edge, and every undirected edge corresponds to a reversible edge.

The \code{dag.m} can be provided using a formula statement (similar to glm). A typical formula is \code{ ~ node1|parent1:parent2 + node2:node3|parent3}. The formula statement have to start with \code{~}. In this example, node1 has two parents (parent1 and parent2). node2 and node3 have the same parent3. The parents names have to exactly match those given in \code{name}. \code{:} is the separator between either children or parents, \code{|} separates children (left side) and parents (right side), \code{+} separates terms, \code{.} replaces all the variables in \code{name}.
}


\value{A matrix giving the PDAG.}

\references{
West, Douglas Brent. Introduction to graph theory. Vol. 2. Upper Saddle River: Prentice hall, (2001).

Further information about \bold{abn} can be found at:\cr
  \url{http://www.r-bayesian-networks.org}}

\author{Gilles Kratzer}


\examples{
dag <- matrix(data = c(0,0,0,
                          1,0,0,
                          1,1,0),nrow = 3,ncol = 3)
dist <- list(a="gaussian",b="gaussian",c="gaussian")     
colnames(dag) <- rownames(dag) <- names(dist)
                          
essentialGraph(dag.m = dag)
}

\keyword{DAG, PDAG}
