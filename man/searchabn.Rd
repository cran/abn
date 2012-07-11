%                               -*- Mode: Rd -*- 
% fitbn.Rd --- 
% Author          : Fraser Lewis
% Created On      : 
% Last Modified By: Fraser Lewis
% Last Modified On: 
% Update Count    : 
% Status          : Unknown, Use with caution!
% 

\name{searchabn}
\alias{searchabn}
\encoding{latin1}
%- Also NEED an `\alias' for EACH other topic documented here.

\title{Search for a locally optimal additive Bayesian network model for continuous and/or binary data}

\description{Implements a local stepwise hill-climbing heuristic searching from an initial random DAG}


\usage{
searchabn (data.df,banned.m,retain.m,start.m,
                     hyper.params=list("mean"=c(0),
                                       "sd"=c(sqrt(1000)),
                                       "shape"=c(0.001),
                                       "scale"=c(1/0.001)),
           max.parents=NULL,init.permuts=0,db.size=10000,max.iters=100,
           epsabs=1e-7,error.verbose=FALSE,enforce.db.size=TRUE,std=TRUE)

}

%- maybe also `usage' for other objects documented here.
\arguments{
  \item{data.df}{a data frame containing the data used for learning each network}
  \item{banned.m}{a matrix defining which arcs are to be banned in the network search, see details}
  \item{retain.m}{a matrix defining which arcs are to be banned in the network search, see details}
  \item{start.m}{a matrix defining which arcs are to be banned in the network search, see details} 
  \item{hyper.params}{prior means and standard deviations in Gaussian prior distribution for regression coefficients, see details}
  \item{max.parents}{maximum number of parents allowed per node in the network search, applies to all nodes}
  \item{init.permuts}{number of permutations applied to an initially empty network from which the search heuristic then commences, see details.} 
  \item{db.size}{the number of individual node score to store so they do not need to be re-calculated during network score calculations.}
  \item{max.iters}{maximum number of iterations used to find root in Laplace approximation}
  \item{epsabs}{minimum absolute error of root estimate in Laplace approximation}
  \item{error.verbose}{logical, if TRUE then extra information is printed to the terminal about any numerical convergence issues in the Laplace approximation}
  \item{enforce.db.size}{logical, if TRUE then if the db.size is exceeded then the search is killed off}
  \item{std}{logical, if TRUE then each variable which is not a factor is standardised to zero mean and standard deviation one} 
}

\details{
 The procedure starts with the DAG given in start.m and randomly adds arcs to this, the number of which is set by init.permuts subject to i) that a cycle is not created; ii) that no node has more than max.parents parents; and iii) the arc is not prohibited by being in banned.m. The choice of K2 or BDeu metric is explained in \code{fitbn}. A each stage in the search the current DAG is permuted by i) adding an arc; ii) removing an arc; and iii) reversing an arc; where these permutations are all subject to the previous three constraints in addition that no arcs in retain.m can be removed. The search progresses to the permuted network with the maximum score (provided that it is greater than that of the current DAG, else the search terminates). Note that only minimal checks are done on start.m and retain.m to check for inconsistencies and retain.m should be included in start.m. See \code{fitbn} and dag.m for the format of the matrices start.m, banned.m and retain.m.  The hyper.parameters are explained in \code{fitabn}.
}


\value{A list where the first component is called scores and is a vector of the network scores for each of the networks found at each iteration of the search. Each remaining entry in the list is a matrix corresponding to the DAG at each stage of the search, the ordering of the scores matches the ordering of the DAGs}


\references{Heckerman, D. and Geiger, D. and Chickering, D. M. (1995). Learning Bayesian Networks - The Combination of Knowledge and Statistical-Data, Machine Learning, 20 (3), 197-243.

Further information about \bold{abn} can be found at:\cr
  \url{http://www.vetepi.uzh.ch/research/bgm.html}
}

\author{
 Fraser Lewis \email{fraseriain.lewis@uzh.ch}
}

\seealso{\code{\link{fitbn}}, \code{\link{hillsearchbn}}}

\examples{
cts.nodes<-c(2,5,7,8,13,14,16,17,22,23,24,25,29,30,31,33);
var33.cts<-var33[,cts.nodes];#drop categorical nodes
mydag<-matrix(c(
                 0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0, #v2
                 0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0, #v5  
                 0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0, #v7  
                 0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0, #v8 
                 0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0, #v13  
                 0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0, #v14  
                 0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0, #v16  
                 0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0, #v17  
                 0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0, #v22 
                 0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0, #v23
                 0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0, #v24 
                 0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0, #v25 
                 0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0, #v29 
                 0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0, #v30 
                 0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0, #v31 
                 0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0  #v33 
              ),byrow=TRUE,ncol=16); 
colnames(mydag)<-rownames(mydag)<-names(var33.cts);#set names
banned.cts<-matrix(rep(0,dim(var33.cts)[2]^2),ncol=dim(var33.cts)[2]);
colnames(banned.cts)<-rownames(banned.cts)<-names(var33.cts);#set names
retain.cts<-matrix(rep(0,dim(var33.cts)[2]^2),ncol=dim(var33.cts)[2]);
colnames(retain.cts)<-rownames(retain.cts)<-names(var33.cts);#set names
start.cts<-matrix(rep(0,dim(var33.cts)[2]^2),ncol=dim(var33.cts)[2]);
colnames(start.cts)<-rownames(start.cts)<-names(var33.cts);#set names
#
myres.add<-searchabn(data.df=var33.cts,
                 banned.m=banned.cts,
                 retain.m=retain.cts,
                 start.m=start.cts,
                 max.parents=2, 
                 init.permuts=0,db.size=1000,error.verbose=TRUE);
}

\keyword{models}
