% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/fitabn_bayes.R
\name{forLoopContentFitBayes}
\alias{forLoopContentFitBayes}
\title{Regress each node on its parents.#'}
\usage{
forLoopContentFitBayes(
  child = NULL,
  dag = NULL,
  data.df = NULL,
  var.types = NULL,
  grouped.vars = NULL,
  group.ids = NULL,
  control = NULL,
  INLA.marginals = NULL,
  verbose = NULL,
  force.method = NULL,
  data.dists = NULL,
  mymodes = NULL,
  error.code = NULL,
  hessian.accuracy = NULL,
  mymargs = NULL
)
}
\arguments{
\item{child}{integer of node to be regressed}

\item{var.types}{vector of numeric encoding of distribution types. See \code{get.var.types(data.dists)}}

\item{INLA.marginals}{vector of logicals indicating which nodes are to be fitted using INLA}

\item{mymodes}{Empty list of modes for each node}

\item{error.code}{Empty element of error codes for each node}

\item{hessian.accuracy}{Empty element of hessian accuracies for each node}

\item{mymargs}{Empty list of marginals for each node}
}
\value{
list of mlik, modes, marginals, error codes, hessian accuracies and a logical if INLA was used for each node.
}
\description{
Regress each node on its parents.#'
}
\keyword{internal}
