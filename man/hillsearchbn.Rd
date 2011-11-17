%                               -*- Mode: Rd -*- 
% fitbn.Rd --- 
% Author          : Fraser Lewis
% Created On      : 
% Last Modified By: Fraser Lewis
% Last Modified On: 
% Update Count    : 
% Status          : Unknown, Use with caution!
% 

\name{hillsearchbn}
\alias{hillsearchbn}
\encoding{latin1}
%- Also NEED an `\alias' for EACH other topic documented here.

\title{Repeated searches for locally optimal Bayesian Networks}

\description{Implements repeated local stepwise hill-climbing heuristic}

\usage{
hillsearchbn (data.df,banned.m,retain.m,start.m,prior.obs.per.node=NULL,
              useK2=FALSE,max.parents=NULL,init.permuts=0,num.searches=1,
              db.size=10000,localdb=TRUE,timing=TRUE,enforce.db.size=TRUE)

}

%- maybe also `usage' for other objects documented here.
\arguments{
  \item{data.df}{a data frame containing the data used for learning the network}
  \item{banned.m}{a matrix defining which arcs are to be banned in the network 
                  search, see details}
  \item{retain.m}{a matrix defining which arcs are to be banned in the network 
                  search, see details} 
  \item{start.m}{a matrix defining which arcs are to be banned in the network search, see details}               
  \item{prior.obs.per.node}{imaginary database size, see details}
  \item{useK2}{logical, choose either K2 metric or BDeu metric, if FALSE BDeu metric
               is used, see details}
  \item{max.parents}{maximum number of parents allowed per node in the network search,
                     applies to all nodes}
  \item{init.permuts}{number of permutations applied to an initially empty network from
                      which the search heuristic then commences, see details.}
  \item{num.searches}{number of times to run the individual local search}
  \item{db.size}{the number of individual node score to store so they do not need to be re-calculated during network score calculations.}
\item{localdb}{should the node cache be reset after each individual search - optimal choice is data dependent}
\item{timing}{output elapsed CPU time taken during each individual search}
\item{enforce.db.size}{logical, if TRUE then if the db.size is exceeded then the search is killed off}
}

\details{
  This runs same search as in \code{searchbn} but repeated num.searches times to produce
  a collection of locally optimal networks. See vignette for more details.

}


\value{A list with entries: init.score; final.score; init.dag; final.dag. The first two are vectors
       containing the scores of the network from which each search commenced and the score of the
       best network found. The initial networks (DAGs) in each of the num.searches are stored as 
       a list of matrices in init.dag, and similarly all the best networks (DAGs) found are stored
       in a list of matrices in final.dag.}


\references{Heckerman, D. and Geiger, D. and Chickering, D. M. (1995). Learning Bayesian Networks - The Combination of Knowledge and Statistical-Data, Machine Learning, 20 (3), 197-243.

Lewis, F. I., Brulisauer, F. and Gunn, G. J. (2011). Structure discovery in Bayesian networks: An analytical tool for analysing complex animal health data Preventive Veterinary Medicine,  100, 109-115.

Further information about \bold{abn} can be found at:\cr
  \url{http://www.vetepi.uzh.ch/}. 
}

\author{
  Fraser Lewis \email{fraseriain.lewis@uzh.ch}
}

\seealso{\code{\link{fitbn}}, \code{\link{searchbn}}}


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
set.seed(10000);## only affects init.permuts
start.list<-list();
n.searches<-2;#example only - must be *much* larger in practice
for(i in 1:n.searches){start.list[[i]]<-retain.cat;} ## empty networks
myres<-hillsearchbn(data.df=var33.cat,
                banned.m=banned.cat,
                retain.m=retain.cat,
                start.m=start.list,
                useK2=TRUE,max.parents=2,init.permuts=0,num.searches=2);

}

\keyword{models}
