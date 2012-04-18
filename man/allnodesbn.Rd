%                               -*- Mode: Rd -*- 
% fitbn.Rd --- 
% Author          : Fraser Lewis
% Created On      : 
% Last Modified By: Fraser Lewis
% Last Modified On: 
% Update Count    : 
% Status          : Unknown, Use with caution!
% 

\name{allnodesbn}
\alias{allnodesbn}
\encoding{latin1}
%- Also NEED an `\alias' for EACH other topic documented here.

\title{Fit all possible parent combinations for Bayesian network}

\description{For use with exact order searches, \code{findmostprobablebn} and \code{getposteriorfeaturesbn}. Calculates network scores for all possible parent combinations for each node subject to the parent limit max.parents and this information is then used in the order based exhaustive search. This is for use with binary or multinomial observations only and fits the usual conjugate BN models as in \code{fitbn}}

\usage{
allnodesbn (data.df,prior.obs.per.node=NULL,useK2=FALSE,max.parents=NULL,
           all.nodes=TRUE, which.nodes=NULL, verbose=FALSE)
}
    
%- maybe also `usage' for other objects documented here.
\arguments{
  \item{data.df}{a data frame containing the data used for learning the network}
  \item{prior.obs.per.node}{imaginary database size, see details}
  \item{useK2}{logical, choose either K2 metric or BDeu metric, if FALSE BDeu metric is used}
  \item{max.parents}{numeric, the maximum number of parents allowed per node}
  \item{all.nodes}{logical, should all parent combinations for all nodes be calculated}
  \item{which.nodes}{numeric vector of node indexes denoting which columns in data.df to include, only used if all.nodes = FALSE}  
  \item{verbose}{logical, if TRUE then parameters estimates are printed for each and every parent combination at each node.} 
  }


\details{
  The procedure fits a Bayesian network model - as per \code{fitbn} - to each node in the data and for all possible parent combinations. This is just a helper function for use with \code{getposteriorfeaturesbn} and \code{findmostprobablebn}. It is designed to allow the calculation across all parents and all nodes to be split into separate computational units, for example to exploit parallel processing. Although, in practice subsequent calls to \code{findmostprobablebn} or \code{getposteriorfeaturesbn} will take vastly more computational time than this function.

}


\value{A list with four elements, a vector of node indexes (from 1 not zero), a matrix where each row is a parent combination for the node in node indexes, score for each parent combination, a flag as to whether a parent limit has been used}

\references{

  Further information about \bold{abn} can be found at:\cr
  \url{http://www.vetepi.uzh.ch/research/bgm.html}
}

\author{
  Fraser Lewis \email{fraseriain.lewis@uzh.ch}
}


\examples{

## See examples in ?getposteriorfeaturesbn


}

\keyword{models}
