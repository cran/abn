%                               -*- Mode: Rd -*- 
% fitbn.Rd --- 
% Author          : Fraser Lewis
% Created On      : 
% Last Modified By: Fraser Lewis
% Last Modified On: 
% Update Count    : 
% Status          : Unknown, Use with caution!
% 

\name{getposteriorfeaturesbn}
\alias{getposteriorfeaturesbn}
\encoding{latin1}
%- Also NEED an `\alias' for EACH other topic documented here.

\title{Find the posterior probability of structural features}

\description{Find the posterior probability for specific structural features, such as arcs, or sub-networks of arcs, in a Bayesian network using order based estimation. }

\usage{
getposteriorfeaturesbn(local.scores,feature="arc",child=NULL, 
                       parent=NULL, featuredefn=NULL,
                       offset=1)

}
    
%- maybe also `usage' for other objects documented here.
\arguments{
  \item{local.scores}{output from allnodesbn or allnodesabn}
  \item{feature}{either "arc", "all" or "custom", see details}
  \item{child}{index of child node - the column index (counting from 1) in the data used in local.scores.}
  \item{parent}{index of parent node - the column index (counting from 1) in the data used in local.scores}
  \item{featuredefn}{a vector of ones or zeros the same length as the first vector in local.score, see details}
  \item{offset}{numeric, a multiplication factor which may help avoid numeric underflow errors when exponentiating logarithms of probabilities}
  }


\details{
  The procedure runs the exact order based structure discovery approach of Koivisto (2004). The local.score is the node cache, and feature then provides information on how to process this cache, e.g. what features to average over. To get the posterior probability for a feature first run this function for the feature, for example an arc. This function returns an unstandardized ``probability'' which must be then divided by an appropriate denominator, which requires running this function again but using "all" for the feature argument. See the example below. The posterior probability for an arc is the sum of marginal likelihoods for all DAGs (grouped into orders) which contain this arc, divided by the total sum of marginal likelihoods. The featuredefn is a vector which should be the same length as the parent combinations for each node in local.scores and is simply a flag which denotes which local log marginal likelihood to include when considering all possible parent combinations. For example this will be a vector comprising of ones and zeros, where each one corresponds to a parent combination in which the arc (or joint collection of arcs) is present. See the abn vignette for more information on specifying featuredefn.

As this algorithm requires the summing of marginal likelihood values then underflow is a possibility when inverting the log marginal likelihood of each parent combination for each node. In some cases supplying a multiplicative offset may avoid this, for example using exp(-mean(nodescores)) for the nodescore in local.score. If no explicit offset is given then this defaults to one.  

}


\value{A number, the unstandardized posterior probability of the structural feature}

\references{ Koivisto, M. V. (2004). Exact Structure Discovery in Bayesian Networks, Journal of Machine Learning Research, vol 5, 549-573.

  Further information about \bold{abn} can be found at:\cr
  \url{http://www.vetepi.uzh.ch/research/bgm.html}
}

\author{
  Fraser Lewis \email{fraseriain.lewis@uzh.ch}
}


\examples{
pigs<-pigs.1par[1:50,c(2,3)];#take a subset of observations and variables
## step 1. find all the local scores within each node - this creates a "cache" of scores
## n.b. using a parent limit of 1 e.g we only consider DAGs with at most one arc per node
mynodes.add<-allnodesabn(data.df=pigs,max.parents=1); 
## step 2. search across all orders/DAGs for the given structural feature
# specifically calculate the get posterior prob of arc 1->2 
p.2.1<-getposteriorfeaturesbn(local.scores=mynodes.add,feature="arc",child=2,parent=1);
print(p.2.1);
## step 3. calculate the denominator - across all structural features (see text)
denom<-getposteriorfeaturesbn(local.scores=mynodes.add,feature="all");
print(denom);
## we now get the posterior probability for arc 1->2
print(p.2.1/denom);

#for more examples see vignette

}

\keyword{models}
