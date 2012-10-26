% abn-internal.Rd --- 
% Author           : Fraser Lewis
% Last modified on : 03/08/2012 

\name{abninla-internal}

\alias{check.data}
\alias{check.dag}
\alias{calc.node.mlike.inla}

\title{abn internal functions}

\description{
  These are functions for either internal use or error checking and are not meant to be called by the user.
}

\usage{
check.data(data.df=NULL,data.dists=NULL, ntrials.loc=NULL, exposure.loc=NULL, use.inla=NULL,group.var=NULL)
check.dag(dag.m,data.df, is.ban.matrix,use.inla,group.var)
tidy.cache(thecache)
}

\details{
check.data() performs basic error checking on data and function arguments, while check.dag() is similar but for the DAG definition and includes a check for cycles. tidy.cache() checks to see if there are any nodes which could not be fitted and then silently removes then from the cache. 

}

\author{
  Fraser Lewis \email{fraseriain.lewis@uzh.ch}.
}


\keyword{internal}
