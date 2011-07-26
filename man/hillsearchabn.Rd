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

\title{Repeated searches for locally optimal additive Bayesian Networks}

\description{Implements repeated local stepwise hill-climbing heuristic}

\usage{
hillsearchabn (data.df,banned.m,hyper.params=list("mean"=c(0),
              "var"=c(1000)), 
              max.parents=NULL,init.permuts=0,num.searches=1)

}

%- maybe also `usage' for other objects documented here.
\arguments{
  \item{data.df}{a data frame containing the data used for learning the network}
  \item{banned.m}{a matrix defining which arcs are to be banned in the network 
                  search, see details}
 \item{hyper.params}{prior mean and var in Gaussian prior distribution for regression coefficients on logit scale, see details}
  \item{max.parents}{maximum number of parents allowed per node in the network search,
                     applies to all nodes}
  \item{init.permuts}{number of permutations applied to an initially empty network from
                      which the search heuristic then commences, see details.}
  \item{num.searches}{number of times to run the individual local search} 
}

\details{
  This runs same search as in \code{searchabn} but repeated num.searches times to produce
  a collection of locally optimal networks. It is structurally similar to \code{hillsearchbn}.

}


\value{A list with entries: init.score; final.score; init.dag; final.dag. The first two are vectors
       containing the scores of the network from which each search commenced and the score of the
       best network found. The initial networks (DAGs) in each of the num.searches are stored as 
       a list of matrices in init.dag, and similarly all the best networks (DAGs) found are stored
       in a list of matrices in final.dag.}


\references{Heckerman, D. and Geiger, D. and Chickering, D. M. (1995). Learning Bayesian Networks
            - The Combination of Knowledge and Statistical-Data, Machine Learning, 20 (3), 197-243.

 Further information about \bold{abn} can be found at:\cr
  \url{http://www.vetepi.uzh.ch/}
}

\author{
  Fraser Lewis \email{fraseriain.lewis@uzh.ch}
}

\seealso{\code{\link{searchbn}}, \code{\link{searchabn}}}


\examples{
data(sim10var);## 10 categorical variables

##  a matrix of arcs to be banned from the search process 
##  - do not ban any arcs - empty matrix
banned<-matrix(c(      0,0,0,0,0,0,0,0,0,0, ## X1
                       0,0,0,0,0,0,0,0,0,0, ## X2
                       0,0,0,0,0,0,0,0,0,0, ## X3
                       0,0,0,0,0,0,0,0,0,0, ## X4
                       0,0,0,0,0,0,0,0,0,0, ## X5
                       0,0,0,0,0,0,0,0,0,0, ## X6
                       0,0,0,0,0,0,0,0,0,0, ## X7
                       0,0,0,0,0,0,0,0,0,0, ## X8
                       0,0,0,0,0,0,0,0,0,0, ## X9
                       0,0,0,0,0,0,0,0,0,0),## X10
                       byrow=TRUE,ncol=10);
colnames(banned)<-names(sim10varadd);# must set names
rownames(banned)<-names(sim10varadd);# must set names

## run the search 
## only consider networks with 3 or less parents per node
\dontrun{myres<-hillsearchabn(data.df=sim10varadd,banned.m=banned, 
                 hyper.params=list(
                 mean=c(0,0,0,0,0,0,0,0,0,0,0),
                 var=c(100,100,100,100,100,100,100,100,100,100,100)),
                 max.parents=3,init.permuts=0, num.searches=2);}

}

\keyword{models}
