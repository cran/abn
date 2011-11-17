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

## see abn vignette for examples 

}

\keyword{models}
