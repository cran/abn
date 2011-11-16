%                               -*- Mode: Rd -*- 
% fitbn.Rd --- 
% Author          : Fraser Lewis
% Created On      : 
% Last Modified By: Fraser Lewis
% Last Modified On: 
% Update Count    : 
% Status          : Unknown, Use with caution!
% 

\name{searchbn}
\alias{searchbn}
\encoding{latin1}
%- Also NEED an `\alias' for EACH other topic documented here.

\title{Search for a locally optimal Bayesian Network for multinomial data}

\description{Implements a local stepwise hill-climbing heuristic searching from an initial random DAG}


\usage{
searchbn (data.df,banned.m,retain.m,start.m,prior.obs.per.node=NULL,
          useK2=FALSE,max.parents=NULL,init.permuts=0,db.size=10000,
          enforce.db.size=TRUE)

}

%- maybe also `usage' for other objects documented here.
\arguments{
  \item{data.df}{a data frame containing the data used for learning each network}
  \item{banned.m}{a matrix defining which arcs are to be banned in the network search, see details}
  \item{retain.m}{a matrix defining which arcs are to be banned in the network search, see details}
 \item{start.m}{a matrix defining which arcs are to be banned in the network search, see details}
  \item{prior.obs.per.node}{imaginary database size, see details}
  \item{useK2}{logical, choose either K2 metric or BDeu metric, if FALSE BDeu metric is used, see details}
  \item{max.parents}{maximum number of parents allowed per node in the network search, applies to all nodes}
  \item{init.permuts}{number of permutations applied to an initially empty network from which the search heuristic then commences, see details.} 
  \item{db.size}{the number of individual node score to store so they do not need to be re-calculated during network score calculations.} 
 \item{enforce.db.size}{logical, if TRUE then if the db.size is exceeded then the search is killed off}
}

\details{
  The procedure starts with the DAG given in start.m and randomly adds arcs to this, the number of which is set by init.permuts subject to i) that a cycle is not created; ii) that no node has more than max.parents parents; and iii) the arc is not prohibited by being in banned.m. The choice of K2 or BDeu metric is explained in \code{fitbn}. A each stage in the search the current DAG is permuted by i) adding an arc; ii) removing an arc; and iii) reversing an arc; where these permutations are all subject to the previous three constraints in addition that no arcs in retain.m can be removed. The search progresses to the permuted network with the maximum score (provided that it is greater than that of the current DAG, else the search terminates). Note that only minimal checks are done on start.m and retain.m to check for inconsistencies and retain.m should be included in start.m. See \code{fitbn} and dag.m for the format of the matrices start.m, banned.m and retain.m.  

}


\value{A list where the first component is called scores and is a vector of the network scores for each of the networks found at each iteration of the search. Each remaining entry in the list is a matrix corresponding to the DAG at each stage of the search, the order of the scores matches the order of the DAGs}


\references{Heckerman, D. and Geiger, D. and Chickering, D. M. (1995). Learning Bayesian Networks - The Combination of Knowledge and Statistical-Data, Machine Learning, 20 (3), 197-243.

Lewis, F. I., Brulisauer, F. and Gunn, G. J. (2011). Structure discovery in Bayesian networks: An analytical tool for analysing complex animal health data Preventive Veterinary Medicine,  100, 109-115.

Further information about \bold{abn} can be found at:\cr
  \url{http://www.vetepi.uzh.ch/}
}

\author{
 Fraser Lewis \email{fraseriain.lewis@uzh.ch}
}

\seealso{\code{\link{fitbn}}, \code{\link{hillsearchbn}}}

\examples{
bin.nodes<-c(1,3,4,6,9,10,11,12,15,18,19,20,21,26,27,28,32); 
var33.cat<-var33[,bin.nodes];#categorical nodes only
mydag<-matrix(c(
                 0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0, #v1
                 0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0, #v3
                 0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0, #v4  
                 0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0, #v6  
                 0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0, #v9  
                 0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0, #v10  
                 0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0, #v11  
                 0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0, #v12  
                 0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0, #v15  
                 0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0, #v18 
                 0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0, #v19
                 0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0, #v20 
                 0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0, #v21 
                 0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0, #v26 
                 0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0, #v27 
                 0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0, #v28 
                 0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0  #v32 
              ),byrow=TRUE,ncol=17); 
colnames(mydag)<-rownames(mydag)<-names(var33.cat);#set names
## create empty DAGs
banned.cat<-matrix(rep(0,dim(var33.cat)[2]^2),ncol=dim(var33.cat)[2]);
colnames(banned.cat)<-rownames(banned.cat)<-names(var33.cat);#set names
retain.cat<-matrix(rep(0,dim(var33.cat)[2]^2),ncol=dim(var33.cat)[2]);
colnames(retain.cat)<-rownames(retain.cat)<-names(var33.cat);#set names
start.cat<-matrix(rep(0,dim(var33.cat)[2]^2),ncol=dim(var33.cat)[2]);
colnames(start.cat)<-rownames(start.cat)<-names(var33.cat);#set names
myres<-searchbn(data.df=var33.cat,
                banned.m=banned.cat,
                retain.m=retain.cat,
                start.m=start.cat,
                useK2=TRUE,max.parents=2,init.permuts=0,db.size=1000);
}

\keyword{models}
