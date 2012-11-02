%                               -*- Mode: Rd -*- 
% fitbn.Rd --- 
% Author          : Fraser Lewis
% Created On      : 
% Last Modified By: Fraser Lewis
% Last Modified On: 
% Update Count    : 
% Status          : Unknown, Use with caution!
% 

\name{search.hillclimber}
\alias{search.hillclimber}
\encoding{latin1}
%- Also NEED an `\alias' for EACH other topic documented here.

\title{Find high scoring directed acyclic graphs using heuristic search}

\description{Find high scoring network (DAG) structures using a random re-starts greedy hill-climber heuristic search.}

\usage{
       search.hillclimber(score.cache=NULL,num.searches=1,seed=0, verbose=FALSE, 
                          timing.on=TRUE,start.dag=NULL,trace=FALSE,
                          support.threshold=0.5,create.graph=FALSE)

}  
%- maybe also `usage' for other objects documented here.
\arguments{
  \item{score.cache}{output from \code{buildscorecache()} or \code{buildscorecache.inla()}} 
  \item{num.searches}{number of times to run the search}
  \item{seed}{non-negative integer which sets the seed in the GSL random number generator}
  \item{verbose}{extra output}
  \item{timing.on}{extra output in terms of duration computation}
  \item{start.dag}{a DAG given as a matrix, see details for format, which can be used to explicity provide a starting point for the structural search}
  \item{trace}{logical, and plots the majority consensus network based on all previous searches, requires Rgraphviz and Cairo} 
  \item{support.threshold}{the proportion of search results - each locally optimal DAG - in which each arc must appear to be a part of the consensus network} 
  \item{create.graph}{save the final graph as an R graph object}
  
}


\details{
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


\author{
  Fraser Lewis \email{fraseriain.lewis@uzh.ch}
}


\examples{
\dontrun{
## example 1 
## use built-in simulated data set

mydat<-ex1.dag.data;## this data comes with abn see ?ex1.dag.data

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

#use simple banlist with no constraints
ban<-matrix(c(
		    #   1 2 3 4 5 6    
                        0,0,0,0,0,0,0,0,0,0, # 1 
			0,0,0,0,0,0,0,0,0,0, # 2
			0,0,0,0,0,0,0,0,0,0, # 3 
			0,0,0,0,0,0,0,0,0,0, # 4
			0,0,0,0,0,0,0,0,0,0, # 5 
			0,0,0,0,0,0,0,0,0,0, # 6
                        0,0,0,0,0,0,0,0,0,0, # 7 
			0,0,0,0,0,0,0,0,0,0, # 8
			0,0,0,0,0,0,0,0,0,0, # 9 
			0,0,0,0,0,0,0,0,0,0  # 10     
                                           ),byrow=TRUE,ncol=10);

colnames(ban)<-rownames(ban)<-names(mydat); #names must be set

retain<-matrix(c(
		    #   1 2 3 4 5 6    
                        0,0,0,0,0,0,0,0,0,0, # 1 
			0,0,0,0,0,0,0,0,0,0, # 2
			0,0,0,0,0,0,0,0,0,0, # 3 
			0,0,0,0,0,0,0,0,0,0, # 4
			0,0,0,0,0,0,0,0,0,0, # 5 
			0,0,0,0,0,0,0,0,0,0, # 6
                        0,0,0,0,0,0,0,0,0,0, # 7 
			0,0,0,0,0,0,0,0,0,0, # 8
			0,0,0,0,0,0,0,0,0,0, # 9 
			0,0,0,0,0,0,0,0,0,0  # 10     
                                           ),byrow=TRUE,ncol=10);
colnames(retain)<-rownames(retain)<-names(mydat); #names must be set

## not run because may take some minutes for buildscorecache() 
## parent limits
max.par<-list("b1"=4,"p1"=4,"g1"=4,"b2"=4,"p2"=4,"b3"=4,"g2"=4,"b4"=4,"b5"=4,"g3"=4);
## now build cache

mycache<-buildscorecache(data.df=mydat,data.dists=mydists,
                     dag.banned=ban, dag.retained=retain,max.parents=max.par);

# now peform 1000 greedy searches
heur.res<-search.hillclimber(score.cache=mycache,num.searches=1000,seed=0,verbose=FALSE,
                 timing.on=FALSE);

# repeat but this time have the majority consensus network plotted as the searches progress
heur.res2<-search.hillclimber(score.cache=mycache,num.searches=1000,seed=0,verbose=FALSE,
                              trace=TRUE,timing.on=FALSE);

## for publication quality output for the consensus network use graphviz direct 
tographviz(dag.m=heur.res$consensus,data.df=mydat,data.dists=mydists,outfile="graphcon.dot");
## and then process using graphviz tools e.g. on linux
system("dot -Tpdf -o graphcon.pdf graphcon.dot");
system("evince graphcon.pdf");
## note the .dot file created can be easily edited manually to provide custom shapes, colours etc.


## example 2 - glmm example - but no difference here as the format of the score cache is identical


mydat<-ex3.dag.data[,c(1:5,14)];## this data comes with abn see ?ex1.dag.data

mydists<-list(b1="binomial",
              b2="binomial",
              b3="binomial",
              b4="binomial",
              b5="binomial"
             );
max.par<-3;

mycache.mixed<-buildscorecache(data.df=mydat,data.dists=mydists,group.var="group",
                               cor.vars=c("b1","b2","b3","b4","b5"),
                               max.parents=max.par, which.nodes=c(1:5));

# now peform 1000 greedy searches
heur.res<-search.hillclimber(score.cache=mycache.mixed,num.searches=1000,seed=0,verbose=FALSE,
                             timing.on=FALSE);

# repeat but this time have the majority consensus network plotted as the searches progress
heur.res<-search.hillclimber(score.cache=mycache.mixed,num.searches=1000,seed=0,verbose=FALSE,
                             trace=TRUE,timing.on=FALSE);

## for publication quality out for the consensus network may be better to use graphviz direct 
tographviz(dag.m=heur.res$consensus,data.df=mydat,data.dists=mydists,group.var="group",
           outfile="graphcon.dot");
## and then process using graphviz tools e.g. on linux
system("dot -Tpdf -o graphcon.pdf graphcon.dot");
system("evince graphcon.pdf");

}}

\keyword{models}
