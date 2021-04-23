%                               -*- Mode: Rd -*- 
% plot-abn.Rd --- 
% Author          : Gilles Kratzer
% Created On      : 16.07.2018
% Last Modified By: 
% Last Modified On: 
% Update Count    : 
% Status          : Unknown, Use with caution!
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

\name{compareDag}
\alias{compareDag}
\encoding{latin1}
%- Also NEED an `\alias' for EACH other topic documented here.

\title{Compare two DAGs}

\description{Function that returns multiple graph metrics to compare two DAGs, known as confusion matrix or error matrix.}

\usage{
compareDag(ref, test, node.names = NULL)
}
    
%- maybe also `usage' for other objects documented here.
\arguments{
  \item{ref}{a matrix or a formula statement (see details for format) defining the reference network structure, a directed acyclic graph (DAG). Note that row names must be set or given in \code{node.names} if the DAG is given via a formula statement.}
  \item{test}{a matrix or a formula statement (see details for format) defining the test network structure, a directed acyclic graph (DAG). Note that row names must be set or given in \code{node.names} if the DAG is given via a formula statement.}
  \item{node.names}{a vector of names if the DAGs are given via formula, see details.}
}


\details{

This R function returns standard Directed Acyclic Graph comparison metrics. In statistical classification, those metrics are known as a confusion matrix or error matrix. Those metrics allows visualization of the difference between different DAGs. In the case where comparing TRUTH to learned structure or two learned structures, those metrics allow the user to estimate the performance of the learning algorithm. In order to compute the metrics, a contingency table is computed of a pondered difference of the adjacency matrices od the two graphs. 

The returns metrics are:
TP = True Positive
TN = True Negative
FP = False Positive
FN = False Negative
CP = Condition Positive (ref)
CN = Condition Negative (ref)
PCP = Predicted Condition Positive (test)
PCN = Predicted Condition Negative (test)

True Positive Rate \deqn{=\frac{\sum TP}{\sum CP}}{=\frac{\sum TP}{\sum CP}}
False Positive Rate \deqn{=\frac{\sum FP}{\sum CN}}{=\frac{\sum FP}{\sum CN}}
Accuracy \deqn{=\frac{\sum TP + \sum TN}{Total population}}{=\frac{\sum TP + \sum TN}{Total population}}
G-measure \deqn{\sqrt {{\frac {TP}{TP+FP}}\cdot {\frac {TP}{TP+FN}}}}{\sqrt {{\frac {TP}{TP+FP}}\cdot {\frac {TP}{TP+FN}}}}
F1-Score \deqn{\frac{2 \sum TP}{2 \sum TP + \sum FN + \sum FP}}{\frac{2 \sum TP}{2 \sum TP + \sum FN + \sum FP}}
Positive Predictive Value \deqn{\frac{\sum TP}{\sum PCP}}{\frac{\sum TP}{\sum PCP}}
False Ommision Rate \deqn{\frac{\sum FN}{\sum PCN}}{\frac{\sum FN}{\sum PCN}}
Hamming-Distance: Number of changes needed to match the matrices.

The \code{ref} or \code{test} can be provided using a formula statement (similar to GLM input). A typical formula is \code{ ~ node1|parent1:parent2 + node2:node3|parent3}. The formula statement have to start with \code{~}. In this example, node1 has two parents (parent1 and parent2). node2 and node3 have the same parent3. The parents names have to exactly match those given in \code{node.names}. \code{:} is the separtor between either children or parents, \code{|} separates children (left side) and parents (right side), \code{+} separates terms, \code{.} replaces all the variables in \code{node.names}.
}


\value{A list giving DAGs comparison metrics. The metrics are: True Positive Rate, False Positive Rate, Accuracy, G-measure, F1-Score, Positive Predictive Value, False Omission Rate, and the Hamming-Distance.}


\references{
Sammut, Claude, and Geoffrey I. Webb. (2017). Encyclopedia of machine learning and data mining. Springer.

Further information about \pkg{abn} can be found at:\cr
  \url{http://r-bayesian-networks.org}}

\author{Gilles Kratzer}

\examples{
test.m <- matrix(data = c(0,1,0,
                          0,0,0,
                          1,0,0), nrow = 3, ncol = 3)
    
ref.m <- matrix(data = c(0,0,0,
                          1,0,0,
                          1,0,0), nrow = 3, ncol = 3)
                          
colnames(test.m) <- rownames(test.m) <- colnames(ref.m) <- colnames(ref.m) <- c("a", "b", "c")
                          
compareDag(ref = ref.m, test = test.m)
}
\keyword{utilities}
\concept{DAG}
