%                               -*- Mode: Rd -*- 
% search.heuristic.Rd --- 
% Author          : Gilles Kratzer
% Created On      : 04/07/2018
% Last Modified By: 17/07/2018
% Last Modified On: 
% Update Count    : 
% Status          : Unknown, Use with caution!
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

\name{search.heuristic}
\alias{search.heuristic}
\encoding{latin1}
%- Also NEED an `\alias' for EACH other topic documented here.

\title{A familly of heuristic algorithms that aims at finding high scoring directed acyclic graphs}

\description{A flexible implementation of multiple greedy search algorithms to find high scoring network (DAG)}

\usage{
search.heuristic(score.cache = NULL,
                        score = "mlik",
                        data.dists = NULL,
                        max.parents = 1,
                        num.searches = 1,
               max.steps = 100,
               seed = 42,
               verbose = FALSE,
               start.dag = NULL,
               dag.retained = NULL,
               dag.banned = NULL,
               algo = "hc",
               tabu.memory = 10,
               temperature = 0.9)
}  
%- maybe also `usage' for other objects documented here.
\arguments{
  \item{score.cache}{output from \code{buildscorecache()}.}
  \item{score}{which score should be used to score the network. Possible choices \code{aic, bic, mdl, mlik}.}
  \item{data.dists}{a named list giving the distribution for each node in the network.} 
  \item{max.parents}{a constant giving the maximum number of parents allowed, should be smaller or equal to the constant used with \code{buildscorecache()}.}
  \item{num.searches}{a positive integer giving the number of different search to run, see details.}
  \item{max.steps}{a constant giving the number of search steps per search, see details.}
  \item{seed}{a non-negative integer which sets the seed.}
  \item{verbose}{extra output, see output for details.}
  \item{start.dag}{a DAG given as a matrix, see details for format, which can be used to explicity provide a starting point for the structural search.}
  \item{dag.retained}{a DAG given as a matrix, see details for format. This is necessary if the score.cache was created using an explicit retain matrix, and the same retain matrix should be used here. \code{dag.retained} is used by the algorithm which generates the initial random DAG to ensure that the necessary arcs are retained.} 
  \item{dag.banned}{a DAG given as a matrix, see details for format. This is necessary if the score.cache was created using an explicit banned matrix, and the same ban matrix should be used here. \code{dag.banned} is used by the algorithm which generates the initial random DAG to ensure that the necessary arcs are banned.} 
  \item{algo}{which heuristic algorithm should be used. Possible choices are: \code{hc, tabu, sa}.}
  \item{tabu.memory}{a non-negative integer number to set the memory of the \code{tabu} search.}
  \item{temperature}{a real number giving the update in temperature for the \code{simulated annealing} search algorithm.}
}


\details{

This function contains a flexible implementation of multiple greedy heuristic algorithm particularly well adapted to \code{abn} framework.

This function is 

A flexible implementation of multiple greedy search algorithms to find high scoring network


The hill-climber procedure \code{hc} is similar, but not identical, to the method presented in Heckerman et al. 1995. (Machine Learning, 20, 197-243). Each search begins with a choice of node ordering and then randomly choose a DAG structure or a user given matrix. If a unique matrix is provided no attemp will 

  The procedure runs a greedy hill climbing search similar, but not identical, to the method presented initially in Heckerman et al. 1995. (Machine Learning, 20, 197-243). Each search begins with a randomly chosen DAG structure where this is constructed in such a way as to attempt to choose a DAG uniformly from the vast landscape of possible structures. The algorithm used is as follows: given a node cache (from buildscorecache()), then within this set of all allowed local parent combinations, a random combination is chosen for each node. This is then combined into a full DAG which is then checked for cycles, where this check iterates over the nodes in a random order. If all nodes have at least one child (i.e. at least one cycle is present) then the first node examined has all its children removed, and the check for cycles is then repeated. If this has removed the only cycle present then this DAG is used at the starting point for the search, if not a second node is chosen (randomly) and the process is then repeated until a DAG is obtained. 

  The actual hill climbing algorithm itself differs slightly from that presented in Heckerman et al. as a full cache of all possible local combinations are available. At each hill-climbing step everything in the node cache is considered, in other words all possible single swaps between members of cache currently present in the DAG and those in the full cache. The single swap which provides the greatest increase in goodness of fit is chosen. A single swap here is the removal or addition of any one node-parent combination present in the cache while avoiding a cycle. This means that as well as all single arc changes (addition or removal), multiple arc changes are also considered at each same step, note however that arc reversal in this scheme takes two steps (as this requires first removal of a parent arc from one node and then addition of a parent arc to a different node). The original algorithm perturbed the current DAG by only a single arc at each step but also included arc reversal. The current implementation may not be any more efficient that the original but is arguably more natural given a pre-computed cache of local scores.

  A start DAG may be provided in which case num.searches must equal 1 - this option is really just to provide a local search around a previously identified optimal DAG.

  This function is designed for two different purposes: i) interactive visualisation; and ii) longer batch processing. The former prvoides easy visual "eyeballing" of data in terms of its majority consensus network (or similar threshold), which is a graphical structure which comprises of all arcs which feature in a given proportion (\code{support.threshold}) of locally optimal DAGs already identified during the run. For example, running 1000 searches with trace=TRUE will after each new search plot the current consensus network based on all previous searches. The general hope is that this structure will stabilize - become fixed - relatively quickly, at least for problems with smaller numbers of nodes. Note that libraries Rgraphviz and Cairo are needed when trace=TRUE. When trace=FALSE then there is no graphical output and as such is rather faster. The format of results is identical in each case but note that the choice of random starting networks (the only random part in the algorithm) used when trace=TRUE and trace=FALSE will differ as these use different random number streams (the former uses R's random number generator, whereas the latter uses gsl's).      

}

\value{A list with entires: 
\item{init.score}{a vector giving network score for initial network from which the search commenced}
\item{final.score}{a vector giving the network score for the locally optimal network}
\item{init.dag}{list of matrices, inital DAGs}
\item{final.dag}{list of matrices, locally optimal DAGs}
\item{consensus}{a matrix holding a DAG}
\item{support.threshold}{percentage supported used to create consensus matrix}
}

\references{Lewis FI, McCormick BJJ (2012). Revealing the complexity of health determinants in resource poor settings. \emph{American Journal Of Epidemiology}. DOI:10.1093/aje/KWS183).  

Further information about \bold{abn} can be found at:\cr
  \url{http://www.r-bayesian-networks.org}}


\author{Fraser Ian Lewis}


\examples{
\dontrun{
##############################################
## example: use built-in simulated data set
##############################################

mydat<-ex1.dag.data; ## this data comes with abn see ?ex1.dag.data

## setup distribution list for each node
mydists<-list(b1="binomial",
              p1="poisson",
              g1="gaussian",
              b2="binomial",
              p2="poisson",
              b3="binomial",
              g2="gaussian",
              b4="binomial",
              b5="binomial",
              g3="gaussian"
             );

mycache<-buildscorecache(data.df=mydat,data.dists=mydists
,max.parents=2);

# now peform 100 greedy searches
heur.res<-search.heuristic(score.cache = mycache,score = "bic",data.dists = mydists,
max.parents = 4,start.dag = "random",num.searches = 100,
max.steps = 500,seed = 42,verbose = FALSE,algo = "hc")


plotabn(dag.m = heur.res$dag,data.dists = mydists)
}}

\keyword{models}
