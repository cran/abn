%                               -*- Mode: Rd -*- 
% fitbn.Rd --- 
% Author          : Fraser Lewis
% Created On      : 
% Last Modified By: Fraser Lewis
% Last Modified On: 
% Update Count    : 
% Status          : Unknown, Use with caution!
% 

\name{allnodesabn}
\alias{allnodesabn}
\encoding{latin1}
%- Also NEED an `\alias' for EACH other topic documented here.

\title{Fit all possible parent combinations for additive Bayesian network}

\description{For use with exact order searches, \code{findmostprobablebn} and \code{getposteriorfeaturesbn}. Calculates network scores for all possible parent combinations for each node subject to the parent limit max.parents and this information is then used in the order based exhaustive search. This is for use with binary or mixed binary and Gaussian observations and fits the same form of BN models as in \code{fitabn}}

\usage{
allnodesabn (data.df,hyper.params=list("mean"=c(0),
                                        "sd"=c(sqrt(1000)),
                                        "shape"=c(0.001),
                                        "scale"=c(1/0.001)),
        max.iters=100,epsabs=1e-7, verbose=FALSE,error.verbose=FALSE,
        max.parents=NULL, all.nodes=TRUE, which.nodes=NULL,std=TRUE)
}
    
%- maybe also `usage' for other objects documented here.
\arguments{
  \item{data.df}{a data frame containing the data used for learning the network, binary nodes must be set as factors}
  \item{hyper.params}{prior means and standard deviations in Gaussian prior distribution for regression coefficients, see fitabn}
  \item{max.iters}{maximum number of iterations used to find root in Laplace approximation}
  \item{epsabs}{minimum absolute error of root estimate in Laplace approximation}
  \item{verbose}{logical, if TRUE then extra information is printed to the terminal}
  \item{error.verbose}{logical, if TRUE then extra information is printed to the terminal about any numerical convergence issues in the Laplace approximation}
  \item{max.parents}{numeric, the maximum number of parents allowed per node}
  \item{all.nodes}{logical, should all parent combinations for all nodes be calculated}
  \item{which.nodes}{numeric vector of node indexes denoting which columns in data.df to include, only used if all.nodes = FALSE}
 \item{std}{logical, if TRUE then each variable which is not a factor is standardised to zero mean and standard deviation one}
  }


\details{
  The procedure fits an additive Bayesian network model to each node in the data and for all possible parent combinations. This is just a helper function for use with \code{getposteriorfeaturesbn} and \code{findmostprobablebn}. It is designed to also allow the calculation across all parents and all nodes to be split into separate computational units, for example to exploit parallel processing. Although, in practice subsequent calls to \code{findmostprobablebn} or \code{getposteriorfeaturesbn} will take vastly more computational time than this function.

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
