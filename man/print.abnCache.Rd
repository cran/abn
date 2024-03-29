% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/abn-methods.R
\name{print.abnCache}
\alias{print.abnCache}
\title{Print objects of class \code{abnCache}}
\usage{
\method{print}{abnCache}(x, digits = 3, ...)
}
\arguments{
\item{x}{Object of class \code{abnCache}}

\item{digits}{number of digits of the results.}

\item{...}{additional parameters. Not used at the moment.}
}
\value{
summary statement of the class of \code{abnCache}.
}
\description{
Print objects of class \code{abnCache}
}
\examples{
## Subset of the build-in dataset, see  ?ex0.dag.data
mydat <- ex0.dag.data[,c("b1","b2","g1","g2","b3","g3")] ## take a subset of cols

## setup distribution list for each node
mydists <- list(b1="binomial", b2="binomial", g1="gaussian",
                g2="gaussian", b3="binomial", g3="gaussian")

# Structural constraints
# ban arc from b2 to b1
# always retain arc from g2 to g1

## parent limits
max.par <- list("b1"=2, "b2"=2, "g1"=2, "g2"=2, "b3"=2, "g3"=2)

## now build the cache of pre-computed scores accordingly to the structural constraints
if(requireNamespace("INLA", quietly = TRUE)){
  # Run only if INLA is available
res.c <- buildScoreCache(data.df=mydat, data.dists=mydists,
                         dag.banned= ~b1|b2, dag.retained= ~g1|g2, max.parents=max.par)
print(res.c)
}
}
\concept{DAG}
