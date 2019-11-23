% entropyData.Rd ---
% Author           : Gilles Kratzer
% Created on :       17.07.2018
% Last modification :
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
\name{entropyData}
\alias{entropyData}
\title{Computes an Empirical Estimation of the Entropy from a Table of Counts}

\usage{
entropyData(freqs.table)
}
\arguments{
  \item{freqs.table}{a table of counts.}}
  
\description{
This function empirically estimates the Shannon entropy from a table of counts using the observed frequencies.
}

\details{The general concept of entropy is defined for probability distributions. The \code{entropyData} function estimates empirical entropy from data. The probability is estimated from data using frequency tables. Then the estimates are plug-in in the definition of the entropy to return the so-called empirical entropy. A common known problem of empirical entropy is that the estimations are biased due to the sampling noise. This is also known that the bias will decrease as the sample size increases.}

\value{Shannon's entropy estimate on natural logarithm scale.}

\author{Gilles Kratzer}
\seealso{\code{\link{discretization}}}

\references{
Cover, Thomas M, and Joy A Thomas. (2012). "Elements of Information Theory". John Wiley & Sons.
}

\examples{
## Generate random variable
rv <- rnorm(n = 100, mean = 0, sd = 2)
dist <- list("gaussian")
names(dist) <- c("rv")

## Compute the entropy through discretization
entropyData(discretization(data.df = rv, data.dists = dist,
                        discretization.method = "fd", nb.states = FALSE))
}
\keyword{utilities}
\concept{DAG}
