%                               -*- Mode: Rd -*- 
% infoDag.Rd --- 
% Author          : Gilles Kratzer
% Created On      : 17/07/2018
% Last Modified By: GK (rd file + naming + examples)
% Last Modified On: 
% Update Count    : 
% Status          : Unknown, Use with caution!
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

\name{infoDag}
\alias{infoDag}
\encoding{latin1}
%- Also NEED an `\alias' for EACH other topic documented here.

\title{Compute standard information for a DAG.}

\description{This function returns standard metrics for DAG description. A list that contains the number of nodes, the number of arcs, the average Markov blanket size, the neighbourhood average set size, the parent average set size and children average set size.}


\usage{
infoDag(dag.m=NULL, 
        name=NULL)
}
    
%- maybe also `usage' for other objects documented here.
\arguments{
  \item{dag.m}{a matrix or a formula statement (see details for format) defining the network structure, a directed acyclic graph (DAG). Note that row names must be set or given in \code{data.dists}}
  \item{name}{a vector of names if the DAG is given via formula, see details.}
}

\details{

This function returns a named list with the following entries: the number of nodes, the number of arcs, the average Markov blanket size, the neighbourhood average set size, the parent average set size and children average set size.

The \code{dag.m} can be provided using a formula statement (similar to glm). A typical formula is \code{ ~ node1|parent1:parent2 + node2:node3|parent3}. The formula statement have to start with \code{~}. In this example, node1 has two parents (parent1 and parent2). node2 and node3 have the same parent3. The parents names have to exactly match those given in \code{name}. \code{:} is the separator between either children or parents, \code{|} separates children (left side) and parents (right side), \code{+} separates terms, \code{.} replaces all the variables in \code{name}.
}


\value{A named list that contains following entries: the number of nodes, the number of arcs, the average Markov blanket size, the neighbourhood average set size, the parent average set size and children average set size.}

\references{
West, Douglas Brent. Introduction to graph theory. Vol. 2. Upper Saddle River: Prentice hall, (2001).

Further information about \bold{abn} can be found at:\cr
  \url{http://www.r-bayesian-networks.org}}

\author{Gilles Kratzer}


\examples{

##creating a dag:
    dag <- matrix(data = c(1,0,0,0,1,0,1,0,1), nrow = 3,ncol = 3)
    dist <- list(a="gaussian",b="gaussian",c="gaussian")
    colnames(dag) <- rownames(dag) <- names(dist)
    
    infoDag(dag.m = dag)
}

\keyword{DAG, information metric}
