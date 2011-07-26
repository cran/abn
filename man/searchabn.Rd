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

\title{Search for a locally optimal additive Bayesian Network}

\description{Implements a local stepwise hill-climbing heuristic searching from an initial random DAG}


\usage{
searchabn (data.df,banned.m,
           hyper.params=list("mean"=c(0),"var"=c(1000)),
           max.parents=NULL,init.permuts=0)

}

%- maybe also `usage' for other objects documented here.
\arguments{
  \item{data.df}{a data frame containing the data used for learning each network}
  \item{banned.m}{a matrix defining which arcs are to be banned in the network search, see details}
  \item{hyper.params}{prior mean and var in Gaussian prior distribution for regression coefficients on logit scale, see details}
  \item{max.parents}{maximum number of parents allowed per node in the network search, applies to all nodes}
  \item{init.permuts}{number of permutations applied to an initially empty network from which the search heuristic then commences, see details.} 
}

\details{
  The procedure starts with a randomly constructed DAG created by adding arcs to an initially empty network, the number of arcs added is set by init.permuts subject to i) that a cycle is not created; ii) that no node has more than max.parents parents; and iii) the arc is not prohibited by being in banned.m. The prior is explained in \code{fitabn}. At each stage in the search the current DAG is permuted by i) adding an arc; ii) removing an arc; and iii) reversing an arc; where these permutations are all subject to the previous three constraints. The search progresses to the permuted network with the maximum score (provided that it is greater than that of the current DAG, else the search terminates).  

  In the banned arcs definition, banned.m, each row represents a node in the network, and the columns in each row define the parents for that particular node, see the example below for the specific format.  

}


\value{A list where the first component is called scores and is a vector of the network scores for each of the networks found at each iteration of the search. Each remaining entry in the list is a matrix corresponding to the DAG at each stage of the search, the ordering of the scores matches the ordering of the DAGs}


\references{Heckerman, D. and Geiger, D. and Chickering, D. M. (1995). Learning Bayesian Networks - The Combination of Knowledge and Statistical-Data, Machine Learning, 20 (3), 197-243.

Further information about \bold{abn} can be found at:\cr
  \url{http://www.vetepi.uzh.ch/}
}

\author{
 Fraser Lewis \email{fraseriain.lewis@uzh.ch}
}

\seealso{\code{\link{fitbn}}, \code{\link{hillsearchbn}}}

\examples{
data(sim10varadd);## 10 categorical variables

##  a matrix of arcs to be banned from the search process
# - do not ban any arcs - empty matrix
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

## run the search - init.permuts=0 set the search to start
## from an empty network
## only consider networks with 3 or less parents per node
\dontrun{myres<-searchabn(data.df=sim10varadd,banned.m=banned, 
                 hyper.params=list(
                 mean=c(0,0,0,0,0,0,0,0,0,0,0),
                 var=c(100,100,100,100,100,100,100,100,100,100,100)),
                 max.parents=3,init.permuts=0);}

## get the score of the network from which the search started
\dontrun{print(myres$scores[1]);}## or myres[[1]][1]
## get the DAG of the network from which the search started
\dontrun{print(myres$iter1);}## or myres[[2]]

}

\keyword{models}
