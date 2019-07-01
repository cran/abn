%                               -*- Mode: Rd -*- 
% mostprobable.Rd --- 
% Author          : Fraser Lewis
% Created On      : 
% Last Modified By: Marta Pittavino
% Last Modified On: Gilles Kratzer
% Update Count    : 
% Status          : Unknown, Use with caution!
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

\name{mostprobable}
\alias{mostprobable}
\encoding{latin1}
%- Also NEED an `\alias' for EACH other topic documented here.

\title{Find most probable DAG structure}

\description{Find most probable DAG structure using exact order based approach due to Koivisto and Sood, 2004}

\usage{
mostprobable(score.cache, prior.choice=1, score="mlik", verbose=TRUE) 
}
    
%- maybe also `usage' for other objects documented here.
\arguments{
  \item{score.cache}{object of class \code{abnCache} typically outputted by from \code{buildscorecache()}.} 
  \item{prior.choice}{an integer, 1 or 2, where 1 is a uniform structural prior and 2 uses a weighted prior, see details}
  \item{score}{character giving which network score should be used to select the structure.}
  \item{verbose}{if TRUE then provides some additional output.}
}


\details{
  The procedure runs the exact order based structure discovery approach of Koivisto and Sood (2004) to find the most probable posterior network (DAG). The local.score is the node cache, as created using \code{\link{buildscorecache}} (or manually provided the same format is used). Note that the scope of this search is given by the options used in local.score, for example by restricting the number of parents or the ban or retain constraints given there. 

This routine can take a long time to complete and is highly sensitive to the number of nodes in the network. It is recommended to use this on a reduced data set to get an idea as to the computational practicality of this approach.  In particular, memory usage can quickly increase to beyond what may be available. For additive models, problems comprising up to 20 nodes are feasible on most machines, memory requirements can increase considerably after this, but then so does the run time making this less practical. It is recommended that some form of over-modelling adjustment is performed on this resulting DAG (unless dealing with vast numbers of observations), for example using parametric bootstrapping which is straightforward to implement in MCMC engines such as JAGS or WinBUGS. See the case studies at \url{http://r-bayesian-networks.org} for details.    

The parameter prior.choice determines the prior used within each individual node for a given choice of parent combination. In Koivisto and Sood (2004) p.554 a form of prior is used which assumes that the prior probability for parent combinations comprising of the same number of parents are all equal. Specifically, that the prior probability for parent set G with cardinality |G| is proportional to 1/[n-1 choose |G|] where there are n total nodes. Note that this favours parent combinations with either very low or very high cardinality which may not be appropriate. This prior is used when prior.choice=2. When prior.choice=1 an uninformative prior is used where parent combinations of all cardinalities are equally likely. The latter is equivalent to the structural prior used in the heuristic searches e.g. \code{search.hillclimber}.   

Note that the network score (log marginal likelihood) of the most probable DAG is not returned as it can easily be computed using \code{\link{fitabn}}, see examples below.  

}


\value{An object of class \code{abnMostprobable}, which is a list containing: a matrix giving the DAG definition of the most probable posterior structure, the cache of pre-computed scores and the score used for selection.}

\references{ Koivisto, M. V. (2004). Exact Structure Discovery in Bayesian Networks, Journal of Machine Learning Research, vol 5, 549-573.

Further information about \bold{abn} can be found at:\cr
  \url{http://r-bayesian-networks.org}}

\author{Fraser Iain Lewis}


\examples{
\dontrun{
## Use a subset of a built-in simulated data set
mydat <- ex0.dag.data[,c("b1","b2","g1","g2","p1","p2")]

## Setup distribution list for each node
mydists <- list(b1="binomial", b2="binomial",
              g1="gaussian", g2="gaussian",
              p1="poisson", p2="poisson")

## Parent limits
max.par <- list("b1"=1, "b2"=1, "g1"=1, "g2"=0, "p1"=1, "p2"=2)

## Now build cache with structural constrains:
##   now ban arc from b2 to b1 
##   always retain arc from g2 to g1

mycache <- buildscorecache(data.df=mydat, data.dists=mydists,
                     dag.banned= ~b1|b2, dag.retained= ~g1|g2, max.parents=max.par)

## Now find the globally best DAG
mp.dag <- mostprobable(score.cache=mycache)

## Get the corresponding best goodness-of-fit network score
fitabn(object = mp.dag)$mlik

##############################
## Second example 
##############################

## this data comes with abn see ?ex1.dag.data
mydat <- ex1.dag.data 

## setup distribution list for each node
mydists <- list(b1="binomial", p1="poisson", g1="gaussian", b2="binomial",
              p2="poisson", b3="binomial", g2="gaussian", b4="binomial",
              b5="binomial", g3="gaussian")

## assum no constraints in ban nor retain

## parent limits
max.par <- list("b1"=4,"p1"=4,"g1"=4,"b2"=4,"p2"=4,"b3"=4,"g2"=4,"b4"=4,"b5"=4,"g3"=4)
## now build cache
mycache <- buildscorecache(data.df = mydat, data.dists = mydists, max.parents = max.par)

## Find the globally best DAG
mp.dag <- mostprobable(score.cache=mycache)
fitabn(object=mp.dag)$mlik

## plot the best model
myres <- fitabn(object=mp.dag,create.graph=TRUE)

plotabn(dag.m = mp.dag$dag, data.dists = mydists)#, fitted.values.abn = myres)

## fit the known true DAG
true.dag <- matrix(data=0,ncol=10, nrow = 10)
colnames(true.dag) <- rownames(true.dag) <- names(mydists)

true.dag["p2",c("b1","p1")] <- 1
true.dag["b3",c("b1","g1","b2")] <- 1
true.dag["g2",c("p1","g1","b2")] <- 1
true.dag["b4",c("g1","p2")] <- 1
true.dag["b5",c("g1","g2")] <- 1
true.dag["g3",c("g1","b2")] <- 1

fitabn(dag.m = true.dag, data.df = mydat, data.dists=mydists)$mlik

#################################################################
## example 3 - models with random effects
#################################################################

## this data comes with abn see ?ex3.dag.data
mydat <- ex3.dag.data[,c(1:4,14)] 

mydists <- list(b1="binomial", b2="binomial", b3="binomial",
              b4="binomial")

## This takes a few seconds
mycache.mixed <- buildscorecache(data.df = mydat, data.dists = mydists,
               group.var = "group", cor.vars = c("b1","b2","b3","b4"),
               max.parents=2, which.nodes=c(1:4))

## find the most probable DAG
mp.dag <- mostprobable(score.cache=mycache.mixed)

## get goodness of fit
fitabn(object = mp.dag, data.df = mydat, data.dists = mydists,
      group.var = "group", cor.vars = c("b1","b2","b3","b4"))$mlik
}}

\keyword{models}
\concept{abn}