%                               -*- Mode: Rd -*- 
% search_hillclimber.Rd --- 
% Author          : Fraser Lewis
% Created On      : 
% Last Modified By: Marta Pittavino
% Last Modified On: 06/10/2014 22.07.2018 (typos)
% Update Count    : 
% Status          : Unknown, Use with caution!
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

\name{searchHillclimber}
\alias{search.hillclimber}
\alias{searchHillclimber}
\encoding{latin1}
%- Also NEED an `\alias' for EACH other topic documented here.

\title{Find high scoring directed acyclic graphs using heuristic search.}

\description{Find high scoring network (DAG) structures using a random re-starts greedy hill-climber heuristic search.}

\usage{
searchHillclimber(score.cache = NULL, num.searches = 1, seed = 0,
         verbose = FALSE, timing.on = TRUE, start.dag = NULL,
         support.threshold = 0.5, create.graph = FALSE, dag.retained = NULL)

}  
%- maybe also `usage' for other objects documented here.
\arguments{
  \item{score.cache}{output from \code{buildscorecache()}.} 
  \item{num.searches}{number of times to run the search.}
  \item{seed}{non-negative integer which sets the seed in the GSL random number generator.}
  \item{verbose}{extra output.}
  \item{timing.on}{extra output in terms of duration computation.}
  \item{start.dag}{a DAG given as a matrix, see details for format, which can be used to explicitly provide a starting point for the structural search,}
  \item{support.threshold}{the proportion of search results - each locally optimal DAG - in which each arc must appear to be a part of the consensus network.} 
  \item{create.graph}{save the final graph as an R graph object.}
  \item{dag.retained}{a DAG given as a matrix, see details for format. This is necessary if the score.cache was created using an explicit retain matrix, and the same retain matrix should be used here. dag.retained is used by the algorithm which generates the initial random DAG to ensure that the necessary arcs are retained.} 
  
}


\details{
  The procedure runs a greedy hill climbing search similar, but not identical, to the method presented initially in Heckerman et al. 1995. (Machine Learning, 20, 197-243). Each search begins with a randomly chosen DAG structure where this is constructed in such a way as to attempt to choose a DAG uniformly from the vast landscape of possible structures. The algorithm used is as follows: given a node cache (from buildscorecache()), then within this set of all allowed local parent combinations, a random combination is chosen for each node. This is then combined into a full DAG which is then checked for cycles, where this check iterates over the nodes in a random order. If all nodes have at least one child (i.e. at least one cycle is present) then the first node examined has all its children removed, and the check for cycles is then repeated. If this has removed the only cycle present then this DAG is used at the starting point for the search, if not a second node is chosen (randomly) and the process is then repeated until a DAG is obtained. 

  The actual hill climbing algorithm itself differs slightly from that presented in Heckerman et al. as a full cache of all possible local combinations are available. At each hill-climbing step everything in the node cache is considered, in other words all possible single swaps between members of cache currently present in the DAG and those in the full cache. The single swap which provides the greatest increase in goodness of fit is chosen. A single swap here is the removal or addition of any one node-parent combination present in the cache while avoiding a cycle. This means that as well as all single arc changes (addition or removal), multiple arc changes are also considered at each same step, note however that arc reversal in this scheme takes two steps (as this requires first removal of a parent arc from one node and then addition of a parent arc to a different node). The original algorithm perturbed the current DAG by only a single arc at each step but also included arc reversal. The current implementation may not be any more efficient that the original but is arguably more natural given a pre-computed cache of local scores.

  A start DAG may be provided in which case num.searches must equal 1 - this option is really just to provide a local search around a previously identified optimal DAG.

  This function is designed for two different purposes: i) interactive visualisation; and ii) longer batch processing. The former provides easy visual "eyeballing" of data in terms of its majority consensus network (or similar threshold), which is a graphical structure which comprises of all arcs which feature in a given proportion (\code{support.threshold}) of locally optimal DAGs already identified during the run. The general hope is that this structure will stabilize - become fixed - relatively quickly, at least for problems with smaller numbers of nodes.

}

\value{A list with entries: 
\item{init.score}{a vector giving network score for initial network from which the search commenced}
\item{final.score}{a vector giving the network score for the locally optimal network}
\item{init.dag}{list of matrices, initial DAGs}
\item{final.dag}{list of matrices, locally optimal DAGs}
\item{consensus}{a matrix holding a DAG}
\item{support.threshold}{percentage supported used to create consensus matrix}
}

\references{Lewis, F. I., and McCormick, B. J. J. (2012). Revealing the complexity of health determinants in resource poor settings. \emph{American Journal Of Epidemiology}. DOI:10.1093/aje/KWS183).  

Further information about \pkg{abn} can be found at:\cr
  \url{http://r-bayesian-networks.org}}


\author{Fraser Iain Lewis}

\examples{
\dontrun{
##############################################
## example 1: use built-in simulated data set
##############################################

## this data comes with abn see ?ex1.dag.data
mydat <- ex1.dag.data 

## setup distribution list for each node
mydists <- list(b1="binomial", p1="poisson", g1="gaussian", b2="binomial",
              p2="poisson", b3="binomial", g2="gaussian", b4="binomial",
              b5="binomial", g3="gaussian")

## Build cache may take some minutes for buildscorecache() 
mycache <- buildscorecache(data.df = mydat, data.dists = mydists, 
                           max.parents = 2);

# now peform 100 greedy searches
heur.res <- searchHillclimber(score.cache = mycache,
                 num.searches = 100, timing.on = FALSE)
plotabn(dag.m = heur.res$consensus, data.dists = mydists)

###########################
## example 2 - glmm example
###########################

## this data comes with abn see ?ex1.dag.data
mydat <- ex3.dag.data[,c(1:4,14)] 

mydists <- list(b1="binomial", b2="binomial", b3="binomial",
              b4="binomial")

## This takes a few seconds
mycache.mixed <- buildscorecache(data.df = mydat, data.dists = mydists,
               group.var = "group", cor.vars = c("b1","b2","b3","b4"),
               max.parents=2, which.nodes=c(1:4))

## Now peform 50 greedy searches
heur.res <- searchHillclimber(score.cache = mycache.mixed, num.searches = 50,
                              timing.on = FALSE)
##  Plot the majority consensus network
plotabn(dag.m = heur.res$consensus, data.dists = mydists)

}}

\keyword{models}
\concept{abn}