%                               -*- Mode: Rd -*- 
% fitbn.Rd --- 
% Author          : Fraser Lewis
% Created On      : 
% Last Modified By: Fraser Lewis
% Last Modified On: 
% Update Count    : 
% Status          : Unknown, Use with caution!
% 

\name{hillsearchabn}
\alias{hillsearchabn}
\encoding{latin1}
%- Also NEED an `\alias' for EACH other topic documented here.

\title{Repeated searches for locally optimal additive Bayesian Networks for continuous and/or binary data}

\description{Implements repeated local stepwise hill-climbing heuristic}

\usage{
hillsearchabn (data.df,banned.m,retain.m,start.m,
                 hyper.params=list("mean"=c(0),
                                   "sd"=c(sqrt(1000)),
                                   "shape"=c(0.001),
                                   "scale"=c(1/0.001)),
               max.parents=NULL,init.permuts=0,num.searches=1,db.size=10000,
               localdb=TRUE,timing=TRUE,max.iters=100,epsabs=1e-7,
               error.verbose=FALSE,enforce.db.size=TRUE,std=TRUE)
}

%- maybe also `usage' for other objects documented here.
\arguments{
  \item{data.df}{a data frame containing the data used for learning the network}
  \item{banned.m}{a matrix defining which arcs are to be banned in the network 
                  search, see details}
  \item{retain.m}{a matrix defining which arcs are to be banned in the network 
                  search, see details}
\item{start.m}{a list of a matrices defining starting DAGs, see details}
 \item{hyper.params}{prior means and standard deviations in Gaussian prior distribution for regression coefficients, see details}
  \item{max.parents}{maximum number of parents allowed per node in the network search,
                     applies to all nodes}
  \item{init.permuts}{number of permutations applied to an initially empty network from
                      which the search heuristic then commences, see details.}
  \item{num.searches}{number of times to run the individual local search} 
  \item{db.size}{the number of individual node score to store so they do not need to be re-calculated during network score calculations.}
\item{localdb}{logical, should the node cache be reset after each individual search - optimal choice is data dependent}
\item{timing}{logical, output elapsed CPU time taken during each individual search.}
\item{max.iters}{maximum number of iterations used to find root in Laplace approximation}
 \item{epsabs}{minimum absolute error of root estimate in Laplace approximation}
\item{error.verbose}{logical, if TRUE then extra information is printed to the terminal about any numerical convergence issues in the Laplace approximation}
\item{enforce.db.size}{logical, if TRUE then if the db.size is exceeded then the search is killed off}
\item{std}{logical, if TRUE then each variable which is not a factor is standardised to zero mean and standard deviation one}
 }

\details{
  This runs same search as in \code{searchabn} but repeated num.searches times to produce
  a collection of locally optimal networks. 

}


\value{A list with entries: init.score; final.score; init.dag; final.dag. The first two are vectors
       containing the scores of the network from which each search commenced and the score of the
       best network found. The initial networks (DAGs) in each of the num.searches are stored as 
       a list of matrices in init.dag, and similarly all the best networks (DAGs) found are stored
       in a list of matrices in final.dag.}


\references{Heckerman, D. and Geiger, D. and Chickering, D. M. (1995). Learning Bayesian Networks
            - The Combination of Knowledge and Statistical-Data, Machine Learning, 20 (3), 197-243.

 Further information about \bold{abn} can be found at:\cr
  \url{http://www.vetepi.uzh.ch/research/bgm.html}
}

\author{
  Fraser Lewis \email{fraseriain.lewis@uzh.ch}
}

\seealso{\code{\link{searchbn}}, \code{\link{searchabn}}}


\examples{
mydag<-matrix(c(
0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,#1
0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,#2  
0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,#3  
0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,#4  
0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,#5  
0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,#6  
0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,#7  
0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,#8  
0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,#9  
0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,#10 
0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,#11 
0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,#12 
0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,#13 
0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,#14 
0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,#15 
0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,#16 
0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,#17 
0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,#18 
0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,#19 
0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,#20 
0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,#21 
0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,#22 
0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,#23 
0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,#24 
0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,#25 
0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,#26 
0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,#27 
0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,#28 
0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,#29 
0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,#30
0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,#31
0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,#32 
0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0 #33
                                               ),byrow=TRUE,ncol=33); 
colnames(mydag)<-rownames(mydag)<-names(var33);#set names
## create empty DAGs
banned<-matrix(rep(0,dim(var33)[2]^2),ncol=dim(var33)[2]);
colnames(banned)<-rownames(banned)<-names(var33);#set names
retain<-matrix(rep(0,dim(var33)[2]^2),ncol=dim(var33)[2]);
colnames(retain)<-rownames(retain)<-names(var33);#set names
start<-matrix(rep(0,dim(var33)[2]^2),ncol=dim(var33)[2]);
colnames(start)<-rownames(start)<-names(var33);#set names

set.seed(10000);## only affects init.permuts
start.list<-list();
n.searches<-2;#example only - must be *much* larger in practice
for(i in 1:n.searches){start.list[[i]]<-retain;} ## empty networks
myres<-hillsearchabn(data.df=var33,banned.m=banned,retain.m=retain,
                     start.m=start.list,
                     hyper.params=list(
                     mean=rep(0,dim(var33)[2]+1),
                     sd=rep(sqrt(1000),dim(var33)[2]+1), 
                     shape=rep(0.001,16),## 16 Gaussian nodes
                     scale=rep(1/0.001,16)## 16 Gaussian nodes
                     ),
                     max.parents=2, 
                     num.searches=n.searches,
                     init.permuts=20,db.size=20000,
                     localdb=TRUE,timing=TRUE);


}

\keyword{models}
