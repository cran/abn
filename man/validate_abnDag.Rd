% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/validation.R
\name{validate_abnDag}
\alias{validate_abnDag}
\title{Check for valid DAG of class \code{abnDag}}
\usage{
validate_abnDag(dag, data.df = NULL, returnDag = TRUE, ...)
}
\arguments{
\item{dag}{dag is either a formula, a matrix  or an object of class 'abnDag'}

\item{data.df}{data frame}

\item{returnDag}{if TRUE (default) returns DAG as matrix.}

\item{...}{additional arguments.}
}
\value{
Either TRUE/FALSE or DAG as matrix.
}
\description{
Beside some basic checks, this function also checks for square matrix with no undirected cycles (trivial cycles) and
for no undirected cycles.
}
\details{
Similar to \code{check.valid.dag()}.
}
\keyword{internal}
