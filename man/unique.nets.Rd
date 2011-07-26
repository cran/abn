%                               -*- Mode: Rd -*- 
% fitbn.Rd --- 
% Author          : Fraser Lewis
% Created On      : 
% Last Modified By: Fraser Lewis
% Last Modified On: 
% Update Count    : 
% Status          : Unknown, Use with caution!
% 

\name{uniquenets}
\alias{uniquenets}
\encoding{latin1}
%- Also NEED an `\alias' for EACH other topic documented here.

\title{Find all unique DAGs}

\description{Given a list of matrices each defining a DAG, returns a vector of indexes of the unique DAGs}


\usage{
uniquenets(dags.list)

}

%- maybe also `usage' for other objects documented here.
\arguments{
  \item{dags.list}{a list of matrices where each matrix defines a valid DAG}
}

\details{In the network structure definition, dags.list, each row represents a node in the network, and the columns in each row define the parents for that particular node, see the example below for the specific format.  

}


\value{a vector of indexes}


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
## run the search - init.permuts=0 - this is for illustration only 
## as this forces each search to use the same empty network
## and hence the results of each search will be identical
## only consider networks with 3 or less parents per node
\dontrun{K2.search<-hillsearchbn(data.df=sim10var,banned.m=dagbanned.df,
        prior.obs.per.node=NULL,useK2=TRUE,max.parents=3,init.permuts=0,
        num.searches=2);}

## now find out how many unique searches
## - two random starting networks could theoretically be the same
\dontrun{indexes<-uniquenets(K2.search$init.dag);} 
## only keep results from distinct searches
\dontrun{K2.search$init.score<-K2.search$init.score[indexes];}
\dontrun{K2.search$final.score<-K2.search$final.score[indexes];}
\dontrun{K2.search$init.dag<-K2.search$init.dag[indexes];}
\dontrun{K2.search$final.dag<-K2.search$final.dag[indexes];}

}

\keyword{models}
