%                               -*- Mode: Rd -*-
% plot-abn.Rd ---
% Author          : Gilles Kratzer
% Created On      : 22.12.2016
% Last Modified By: GK (rd file + naming + examples)
% Last Modified On: 01/03/2017 22.07.2018 (typos)
% Update Count    :
% Status          : Unknown, Use with caution!
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

\name{plotabn}
\alias{plotAbn}
\encoding{latin1}
%- Also NEED an `\alias' for EACH other topic documented here.

\title{Plot an ABN graphic}

\description{Plot an ABN DAG using formula statement or a matrix in using Rgraphviz through the graphAM class.}

\usage{
plotAbn(dag, data.dists=NULL, markov.blanket.node=NULL, fitted.values=NULL,
        digits=2, edge.strength=NULL, edge.direction="pc", edge.color="black",
        edge.arrowsize=0.5, edge.fontsize=5, node.fontsize=10,
        node.fillcolor=c("lightblue","brown3","chartreuse3"), node.fillcolor.list=NULL,
        node.shape=c("circle","box","ellipse","diamond"), plot=TRUE , ...)

}

%- maybe also `usage' for other objects documented here.
\arguments{
  \item{dag}{a matrix or a formula statement (see details for format) defining the network structure, a Directed Acyclic Graph (DAG). Note that rownames must be set or given in \code{data.dists}.}
  \item{data.dists}{a named list giving the distribution for each node in the network, see details.}
  \item{markov.blanket.node}{name of variables to display its Markov blanket.}
  \item{fitted.values}{modes or coefficents outputted from \code{\link{fitAbn}}.}
  \item{digits}{number of digits to display the \code{fitted.values}.}
  \item{edge.strength}{a named matrix containing evaluations of edge strength which will change the arcs width (could be Mutual information, p-values, number of bootstrap retrieve samples or the outcome of the \code{\link{link.strength}}).}
   \item{edge.direction}{character giving the direction in which arcs should be plotted, \code{pc} (parent to child) or \code{cp} (child to parent) or \code{undirected}.}
  \item{node.fillcolor}{the colour of the node. Second and third element is used for the Markov blanket and node of the Markov blanket.}
  \item{node.shape}{the shape of the nodes according the Gaussian, binomial, Poisson and multinomial distributions.}
  \item{edge.color}{the colour of the edge.}
  \item{edge.arrowsize}{the thickness of the arrows. Not relevant if \code{arc.strength} is provided.}
  \item{node.fontsize}{the font size of the nodes names.}
  \item{edge.fontsize}{the font size of the arcs fitted values.}
  \item{plot}{logical variable, if set to \code{TRUE} then the graph is plotted.}
  \item{node.fillcolor.list}{the list of node that should be coloured.}
  \item{...}{arguments passed to the plotting function.}
}


\details{

By default binomial nodes are squares, multinoial nodes are empty, Gaussian nodes are circles and poison nodes are ellipses.

The \code{dag} can be provided using a formula statement (similar to glm). A typical formula is \code{ ~ node1|parent1:parent2 + node2:node3|parent3}. The formula statement have to start with \code{~}. In this example, node1 has two parents (parent1 and parent2). node2 and node3 have the same parent3. The parents names have to exactly match those given in \code{name}. \code{:} is the separator between either children or parents, \code{|} separates children (left side) and parents (right side), \code{+} separates terms, \code{.} replaces all the variables in \code{name}.
}


\value{A matrix giving the DAG definition of the most probable posterior structure.}

\references{
Further information about \pkg{abn} can be found at:\cr
  \url{http://r-bayesian-networks.org}}
\author{Gilles Kratzer}
\examples{
#Define distribution list
dist <- list(a="gaussian", b="gaussian", c="gaussian", d="gaussian", e="binomial", f="binomial")

#Define a matrix formulation
arc.strength <- matrix(c(0,0.5,0.5,0.7,0.1,0,
                              0,0,0.3,0.1,0,0.8,
                              0,0,0,0.35,0.66,0,
                              0,0,0,0,0.9,0,
                              0,0,0,0,0,0.8,
                              0,0,0,0,0,0),nrow = 6L, ncol = 6L, byrow = TRUE)

#Naming of the matrix
colnames(arc.strength) <- rownames(arc.strength) <- names(dist)

#Plot from a formula
plotAbn(dag = ~a|b:c:e+b|c:d:f+e|f, data.dist = dist)

#Plot form a matrix
plotAbn(dag = arc.strength, data.dist = dist)

#Arc strength
plotAbn(dag = ~a|b:c:e+b|c:d:f+e|f, data.dist = dist, arc.strength = arc.strength)

#Markov blanket
plotAbn(dag = ~a|b:c:e+b|c:d:f+e|f, data.dists = dist, markov.blanket.node = "e")
}
\keyword{models}
\keyword{hplot}
