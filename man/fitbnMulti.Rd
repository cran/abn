%                               -*- Mode: Rd -*- 
% fitbn.Rd --- 
% Author          : Fraser Lewis
% Created On      : 
% Last Modified By: Fraser Lewis
% Last Modified On: 
% Update Count    : 
% Status          : Unknown, Use with caution!
% 

\name{fitbnMulti}
\alias{fitbnMulti}
\encoding{latin1}
%- Also NEED an `\alias' for EACH other topic documented here.

\title{Fit a list of DAGs to data}

\description{Takes a list of DAGs and fits each one in turn to observed data returning the network score}

\usage{
fitbnMulti (data.df, dag.list, prior.obs.per.node=NULL, useK2=FALSE,
            verbose=FALSE)

}

%- maybe also `usage' for other objects documented here.
\arguments{
  \item{data.df}{a data frame containing the data used for learning the network}
  \item{dag.list}{a list of matrices each defining the network structure, a directed acyclic graph}
  \item{prior.obs.per.node}{imaginary database size}
  \item{useK2}{logical, choose either K2 metric or BDeu metric, if FALSE BDeu metric is used}
  \item{verbose}{logical, if TRUE then parameters estimates are provided for each and every parent combination at each node.} 
}

\details{
  This is the same as \code{fitbn}, the only difference being that a list of DAGs is passed rather than a single DAG. This function is for computational and memory efficiency, the same can be achieved through repeated calls to \code{fitbn} but this is much less efficient}


\value{A vector of network scores.}


\author{
  Fraser Lewis \email{fraseriain.lewis@uzh.ch}
}

\seealso{\code{\link{fitbn}}}

\examples{
data(sim10var);## 10 categorical variables

##  a matrix of arcs to be banned from the search process - do not ban any arcs 
##  - empty matrix
dagbanned.df<-matrix(c(0,0,0,0,0,0,0,0,0,0, ## var1
                       0,0,0,0,0,0,0,0,0,0, ## var2
                       0,0,0,0,0,0,0,0,0,0, ## var3
                       0,0,0,0,0,0,0,0,0,0, ## var4
                       0,0,0,0,0,0,0,0,0,0, ## var5
                       0,0,0,0,0,0,0,0,0,0, ## var6
                       0,0,0,0,0,0,0,0,0,0, ## var7
                       0,0,0,0,0,0,0,0,0,0, ## var8
                       0,0,0,0,0,0,0,0,0,0, ## var9
                       0,0,0,0,0,0,0,0,0,0),## var10
                       byrow=TRUE,ncol=10);

colnames(dagbanned.df)<-names(sim10var);# must set names
rownames(dagbanned.df)<-names(sim10var);# must set names
## run the search 
## only consider networks with 3 or less parents per node
## only use first 100 observations for speed
\dontrun{K2.search<-hillsearchbn(data.df=sim10var[1:100,],banned.m=dagbanned.df, 
         prior.obs.per.node=NULL,useK2=TRUE,max.parents=3,init.permuts=4,
         num.searches=200);}
## now prune each of the 200 locally optimal networks (note - in practice must check
## for duplicate searches - unlikely but not impossible)
\dontrun{prunned<-prunenets(K2.search$final.dag,
         round(0.5*length(K2.search$final.dag)));}#50% threshold for pruning
## only retain arc which are present in at least 50% of all locally optimal networks
## we have now pruned all the networks so we need to re-calculate the scores 
## as arcs may have been removed.
\dontrun{prune.res<-fitbnMulti(data.df=sim10var[1:100,],dag.list=prunned$dags,
         prior.obs.per.node=NULL,useK2=TRUE,verbose=FALSE);}


}

\keyword{models}
