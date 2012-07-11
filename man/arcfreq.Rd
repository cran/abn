%                               -*- Mode: Rd -*- 
% fitbn.Rd --- 
% Author          : Fraser Lewis
% Created On      : 
% Last Modified By: Fraser Lewis
% Last Modified On: 
% Update Count    : 
% Status          : Unknown, Use with caution!
% 

\name{arcfreq}
\alias{arcfreq}
\encoding{latin1}
%- Also NEED an `\alias' for EACH other topic documented here.

\title{Summarise many locally optimal networks}

\description{Summarise many locally optimal networks by frequency counts of arcs}


\usage{
arcfreq(dags.list)

}

%- maybe also `usage' for other objects documented here.
\arguments{
  \item{dags.list}{a list of matrices where each matrix defines a valid DAG}
 % \item{threshold.freq}{the number of DAGs each arc must appear in for it to avoid being pruned}
}

\details{Calculates the frequency of every possible arc across the set of given DAGs. For example arcs which appear in at least 50\% of DAGs would then comprise the majority consensus network.} 
% Based on this prune each network back to a given level, e.g. a consensus of 50\% support removes from each locally optimal network each arc which is not present in at least 50\% of the other locally optimal networks.}
  
\value{a matrix showing the number of locally optimal networks which contain each arc.}


\author{
  Fraser Lewis \email{fraseriain.lewis@uzh.ch}
}



\examples{

## see abn vignette for examples 

}

\keyword{models}
