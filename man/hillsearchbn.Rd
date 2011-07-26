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
hillsearchbn (data.df, banned.m, prior.obs.per.node=NULL, useK2=FALSE,
              max.parents=NULL,init.permuts=0,num.searches=1)

}

%- maybe also `usage' for other objects documented here.
\arguments{
  \item{data.df}{a data frame containing the data used for learning the network}
  \item{banned.m}{a matrix defining which arcs are to be banned in the network 
                  search, see details}
  \item{prior.obs.per.node}{imaginary database size, see details}
  \item{useK2}{logical, choose either K2 metric or BDeu metric, if FALSE BDeu metric
               is used, see details}
  \item{max.parents}{maximum number of parents allowed per node in the network search,
                     applies to all nodes}
  \item{init.permuts}{number of permutations applied to an initially empty network from
                      which the search heuristic then commences, see details.}
  \item{num.searches}{number of times to run the individual local search} 
}

\details{
  This runs same search as in \code{searchbn} but repeated num.searches times to produce
  a collection of locally optimal networks.

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
data(sim10var);## 10 categorical variables

##  a matrix of arcs to be banned from the search process 
##  - do not ban any arcs - empty matrix
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
\dontrun{K2.search<-hillsearchbn(data.df=sim10var,banned.m=dagbanned.df, 
         prior.obs.per.node=NULL, useK2=TRUE,max.parents=3,init.permuts=4,
         num.searches=2);}

}

\keyword{models}
