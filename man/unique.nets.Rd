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

## See abn vignette for examples 

}

\keyword{models}
