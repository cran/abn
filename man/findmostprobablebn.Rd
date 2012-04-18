%                               -*- Mode: Rd -*- 
% fitbn.Rd --- 
% Author          : Fraser Lewis
% Created On      : 
% Last Modified By: Fraser Lewis
% Last Modified On: 
% Update Count    : 
% Status          : Unknown, Use with caution!
% 

\name{findmostprobablebn}
\alias{findmostprobablebn}
\encoding{latin1}
%- Also NEED an `\alias' for EACH other topic documented here.

\title{Find most probable structure}

\description{Find most probable structure using exact order based approach due to Koivisto 2004}

\usage{
       findmostprobablebn(local.scores, data.df=NULL) 

}
    
%- maybe also `usage' for other objects documented here.
\arguments{
  \item{local.scores}{output from allnodesbn or allnodesabn} 
  \item{data.df}{a data frame containing the data used for learning the network. Binary variables should be declared as factors}
 % \item{offset}{numeric, a multiplication factor which may help avoid numeric underflow errors}
  }


\details{
  The procedure runs the exact order based structure discovery approach of Koivisto (2004) to find the most probable posterior network. The local.score is the node cache, as created using \code{allnodesbn} or \code{allnodesabn}. Node that the scope of this search is given by the options used in local.score, for example by restricting the number of parents. Unlike in \code{getposteriorfeaturesbn} numerical underflow is not an issue as this routine only uses maximization and requires no summations which means all computation can be kept on a logarithmic scale.

This routine can take a very long time to complete and is highly sensitive to the number of nodes in the network, and also the number of parents allowed. It is recommended to use this on a reduced data set - fewer variables - to get an idea as to the computational practicality of this approach. For example for more than 14 variables, even using a limit of one parent per node may take days rather than hours to complete on a 3Ghz cpu. Examining variable subsets may be all that is feasible for larger scale problems.  

}


\value{A matrix giving the DAG definition of the most probable posterior structure.}

\references{ Koivisto, M. V. (2004). Exact Structure Discovery in Bayesian Networks, Journal of Machine Learning Research, vol 5, 549-573.

  Further information about \bold{abn} can be found at:\cr
  \url{http://www.vetepi.uzh.ch/research/bgm.html}
}

\author{
  Fraser Lewis \email{fraseriain.lewis@uzh.ch}
}


\examples{
\dontrun{
pigs<-pigs.1par[,c(1:8,12,13)];#all 9011 observations but limit to 10 variables
# using all 13 variables in pigs will take several hours of cpu time
mynodes.add<-allnodesabn(data.df=pigs,max.parents=1); 
## now for most probable network of all DAGs where each node has at
## most one arc
map.1par.10var<-findmostprobablebn(local.scores=mynodes.add,data.df=pigs);
}

## for more examples see vignette

}

\keyword{models}
