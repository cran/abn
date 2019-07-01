% mb.Rd ---
% Author           : Gilles Kratzer
% Created on :       17.07.2018
% Last modification : 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
\name{mb}
\alias{mb}
\title{Compute the Markov blanket}

\description{
This function computes the Markov blanket of a set of node given a DAG (Directed Acyclic Graph).
}


\usage{
mb(dag, node, data.dists=NULL)}
\arguments{
  \item{dag}{a matrix or a formula statement (see details for format) defining the network structure, a directed acyclic graph (DAG).}
  \item{node}{a character vector of the nodes for which the Markov Blanket should be returned.} 
  \item{data.dists}{a named list giving the distribution for each node in the network, see details.}
 
}


\details{
This function return the Markov Blanket of a set of nodes given a DAG. 

The \code{dag} can be provided using a formula statement (similar to glm). A typical formula is \code{ ~ node1|parent1:parent2 + node2:node3|parent3}. The formula statement have to start with \code{~}. In this example, node1 has two parents (parent1 and parent2). node2 and node3 have the same parent3. The parents names have to exactly match those given in \code{name}. \code{:} is the separtor between either children or parents, \code{|} separates children (left side) and parents (right side), \code{+} separates terms, \code{.} replaces all the variables in \code{name}.

}

\examples{
## Defining distribution and dag
dist <- list(a="gaussian", b="gaussian", c="gaussian", d="gaussian", e="binomial", 
               f="binomial")
dag <- matrix(c(0,1,1,0,1,0,
    0,0,1,1,0,1, 
    0,0,0,0,0,0, 
    0,0,0,0,0,0, 
    0,0,0,0,0,1, 
    0,0,0,0,0,0), nrow = 6L, ncol = 6L, byrow = TRUE)
    
colnames(dag) <- rownames(dag) <- names(dist)
    
mb(dag, node = "b", data.dists = dist)
mb(dag, node = "e", data.dists = dist)
mb(dag, node = c("b","e"), data.dists = dist)

}
\keyword{utilities}
\concept{DAG}
