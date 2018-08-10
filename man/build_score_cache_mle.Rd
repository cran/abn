% buildscorecache.mle.Rd --- 
% Author           : Gilles Kratzer
% Created on :       18.07.2018 
% Last modification : 
% Last modification : 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
\name{buildscorecache.mle}

\alias{buildscorecache.mle}

\title{Build a cache of goodness of fit metrics based on Information Theoretic for each node in a DAG, possibly subject to user defined restrictions}

\description{
  Iterates over all valid parent combinations - subject to ban, retain and max.parent limits - for each node, and computes a cache of information theoretic scores. This cache can then be used in different DAG structural search algorithms. 
}

\usage{



buildscorecache.mle(data.df = NULL, 
                                data.dists = NULL, 
                                max.parents = NULL,
                                adj.vars = NULL, 
                                cor.vars = NULL, 
                                dag.banned = NULL,
                                dag.retained = NULL,
                                maxit = 100, 
                                tol = 10^-8,
                                centre = TRUE, 
                                dry.run = FALSE)
}
\arguments{
  \item{data.df}{a data frame containing the data used for learning each node, binary variables must be declared as factors.}
  \item{data.dists}{a named list giving the distribution for each node in the network, see details.} 
  \item{max.parents}{a constant or named list giving the maximum number of parents allowed, the list version allows this to vary per node.}
  \item{adj.vars}{a character vector giving the column names in \code{data.df} for which the network score has to be adjusted for, see details.}
  \item{cor.vars}{a character vector giving the column names in \code{data.df} for which adjustment should be used.}
  \item{dag.banned}{a matrix or a formula statement defining which arcs are not permitted - banned - see details for format. Note that colnames and rownames must be set, otherwise same row/column names as \code{data.df} will be assumed. If set as NULL an empty matrix is assumed.}
  \item{dag.retained}{a matrix or a formula statement (see details for format) defining which arcs are must be retained in any model search, see details for format. Note that colnames and rownames must be set, otherwise same row/column names as \code{data.df} will be assumed. If set as NULL an empty matrix is assumed.}
  \item{maxit}{integer given the maximum number of run for estimating network scores using an Iterative Reweighed Least Square algorithm.}
  \item{tol}{real number giving the minimal tolerance expected to terminate the Iterative Reweighed Least Square algorithm to estimate network score.}
  \item{centre}{logical variable, should the observations in each Gaussian node first be standardised to mean zero and standard deviation one, defaults is TRUE.}
    \item{dry.run}{logical variable, if set to TRUE then a list of the child nodes and parent combinations are returned but without estimation the network score.}

}


\value{A named list containing:
\item{children}{a vector of the child node indexes (from 1) corresponding to the columns in data.df (ignoring any grouping variable)}
\item{node.defn}{a matrix giving the parent combination}
\item{mlik}{log marginal likelihood value for each node combination. If the model cannot be fitted then NaN is returned.}
\item{error.code}{NULL (for compatibility purpose)}
\item{error.code.desc}{NULL (for compatibility purpose)}
\item{hessian.accuracy}{NULL (for compatibility purpose)}
\item{data.df}{a version of the original data (for internal use only in other functions such as \code{\link{mostprobable}}).}
\item{aic}{aic value for each node combination. If the model cannot be fitted then NaN is returned.}
\item{bic}{bic value for each node combination. If the model cannot be fitted then NaN is returned.}
\item{mdl}{mdl value for each node combination. If the model cannot be fitted then NaN is returned.}
}

\details{
This function is used to calculate all individual Information-Theoretic node scores. The possible Information-theoretic based network scores computed in \code{buildscorecache.mle} are the maximum likelihood (mlik, called marginal likelihood in this context as it is computed node wise), the Akaike Information Criteria (aic), the Bayesian Information Criteria (bic) and the Minimum distance Length (mdl). The classical definitions of those metrics are given in Kratzer and Furrer (2018). This function computes a cache that can be fed into a model search algorithm. This function is very similar to \code{\link{fitabn.mle}} - see that help page for details of the type of models used and in particular \code{data.dists} specification - but rather than fit a single complete DAG \code{buildscorecache.mle} iterates over all admissible different parent combinations for each node.  There are three ways to customise the parent combinations through giving a matrix which contains arcs which are not allowed (banned), a matrix which contains arcs which must always be included (retained) and also a general complexity limit which restricts the maximum number of arcs allowed to terminate at a node (its number of parents). In these matrices, \code{dag.banned} and \code{dag.retained}, each row represents a node in the network, and the columns in each row define the parents for that particular node, see the example below for the specific format. If these are not supplied they are assumed to be empty matrices, i.e. no arcs banned or retained. Note that only rudimentary consistency checking is done here and some care should be taken not to provide conflicting restrictions in the ban and retain matrices.

The numerical routines used here are optimized and lighted version of those in \code{\link{fitabn.mle}} (but stay essentially the same) and see that help page for further details.
}

\references{Kratzer, G., Furrer, R., 2018. Information-Theoretic Scoring Rules to Learn Additive Bayesian Network Applied to Epidemiology. Preprint; Arxiv: stat.ML/1808.01126.

  Further information about \bold{abn} can be found at:\cr
  \url{http://www.r-bayesian-networks.org}}

\author{Gilles Kratzer}

\examples{
mydat <- ex0.dag.data[,c("b1","b2","g1","g2","b3","g3")] ## take a subset of cols

## setup distribution list for each node
mydists <- list(b1="binomial",
              b2="binomial",
              g1="gaussian",
              g2="gaussian",
              b3="binomial",
              g3="gaussian"
             )
             
## now build cache of scores (goodness of fits for each node)

res.mle <- buildscorecache.mle(data.df=mydat,data.dists=mydists,max.parents=3)
res.abn <- buildscorecache(data.df=mydat,data.dists=mydists,max.parents=3)

#plot(-res.mle$bic,res.abn$mlik)
}

\keyword{glm, network score}
