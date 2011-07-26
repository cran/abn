%                               -*- Mode: Rd -*- 
% fitbn.Rd --- 
% Author          : Fraser Lewis
% Created On      : 
% Last Modified By: Fraser Lewis
% Last Modified On: 
% Update Count    : 
% Status          : Unknown, Use with caution!
% 

\name{prunenets}
\alias{prunenets}
\encoding{latin1}
%- Also NEED an `\alias' for EACH other topic documented here.

\title{Prune arcs from DAGs according to a threshold criteria}

\description{Summarise many locally optimal networks by pruning}


\usage{
prunenets(dags.list,threshold.freq)

}

%- maybe also `usage' for other objects documented here.
\arguments{
  \item{dags.list}{a list of matrices where each matrix defines a valid DAG}
  \item{threshold.freq}{the number of DAGs each arc must appear in for it to avoid being pruned}
}

\details{Calculate the support across all optimal DAGs found for each arc present in this set. Based on this prune each network back to a given level, e.g. a consensus of 50\% support removes from each locally optimal network each arc which is not present in at least 50\% of the other locally optimal networks.}
  


\value{a list, the first entry, arcs.sum is a matrix showing the number of locally optimal networks which contain each arc. The second entry is called dags, and is a list containing each locally optimal DAG now pruned to the given level.}


\author{
  Fraser Lewis \email{fraseriain.lewis@uzh.ch}
}



\examples{
data(sim10var);## 10 categorical variables

##  a matrix of arcs to be banned from the search process
## - do not ban any arcs - empty matrix
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
## - only use first 100 observations for speed
\dontrun{K2.search<-hillsearchbn(data.df=sim10var[1:100,],
         banned.m=dagbanned.df,
         prior.obs.per.node=NULL,useK2=TRUE,max.parents=3,init.permuts=4,
         num.searches=200);}
## now prune each of the 200 locally optimal networks (note - note checked
##  for unique searches!)
\dontrun{prunned<-prunenets(K2.search$final.dag,
         round(0.5*length(K2.search$final.dag)));}#50% threshold for pruning
\dontrun{print(prunned$arcs.sum);}## the total support for each arc 
##in all the locally optimal networks found
\dontrun{print(prunned$dags[[1]]);}
#the first locally optimal network after pruning
}

\keyword{models}
