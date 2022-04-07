% buildScoreCache.Rd ---
% Author           : Fraser Lewis & Gilles Kratzer
% Created on :       03/08/2012
% Last modification : 01/12/2016
% Last modification : 01/03/2017
% Last modification : 21/05/2019 (GK)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

\name{buildScoreCache}

\alias{buildScoreCache}

\title{Build a cache of goodness of fit metrics for each node in a DAG, possibly subject to user-defined restrictions}

\description{
  Iterates over all valid parent combinations - subject to ban, retain, and \code{max.parent} limits - for each node, or a subset of nodes, and computes a cache of log marginal likelihoods. This cache can then be used in different DAG structural search algorithms.
}

\usage{
buildScoreCache(data.df = NULL, data.dists = NULL, method = "bayes",
   group.var = NULL, adj.vars = NULL, cor.vars = NULL, dag.banned = NULL,
   dag.retained = NULL, max.parents = NULL, which.nodes=NULL,
   defn.res = NULL, centre = TRUE, dry.run = FALSE,
   control = NULL, verbose = FALSE, ...)
}
\arguments{
  \item{data.df}{a data frame containing the data used for learning each node, binary variables must be declared as factors.}
  \item{data.dists}{a named list giving the distribution for each node in the network, see \sQuote{Details}.}
  \item{method}{should a "Bayes" or "mle" approach be used, see \sQuote{Details}.}
  \item{group.var}{only applicable for nodes to be fitted as a mixed
    model (Bayesian) and gives the column name in \code{data.df} of the grouping variable which must be a factor denoting group membership.}
  \item{adj.vars}{a character vector giving the column names in \code{data.df} for which the network score has to be adjusted for, see \sQuote{Details}.}
  \item{cor.vars}{a character vector giving the column names in \code{data.df} for which a mixed model should be used to adjust for within group correlation or pure adjustment.}
  \item{dag.banned}{a matrix or a formula statement (see \sQuote{Details} for format) defining which arcs are not permitted - banned - see \sQuote{Details} for format. Note that colnames and rownames must be set, otherwise same row/column names as data.df will be assumed. If set as NULL an empty matrix is assumed.}
  \item{dag.retained}{a matrix or a formula statement (see \sQuote{Details} for format) defining which arcs are must be retained in any model search, see \sQuote{Details} for format. Note that colnames and rownames must be set, otherwise same row/column names as data.df will be assumed. If set as NULL an empty matrix is assumed.}
  \item{max.parents}{a constant or named list giving the maximum number of parents allowed, the list version allows this to vary per node.}
  \item{which.nodes}{a vector giving the column indices of the variables to be included, if ignored all variables are included.}
  \item{defn.res}{an optional user-supplied list of child and parent combinations, see \sQuote{Details}.}
    \item{centre}{should the observations in each Gaussian node first be standardized to mean zero and standard deviation one, defaults to TRUE.}
  \item{dry.run}{if TRUE then a list of the child nodes and parent combinations are returned but without estimation of node scores (log marginal likelihoods).}
  \item{control}{a list of control parameters. See \code{\link{build.control}} for the names of the settable control values and their effect.}
    \item{verbose}{if TRUE then provides some additional output.}
    \item{\dots}{additional arguments passed for optimization.}
}

\value{A named list of class \code{abnCache}.
\item{children}{a vector of the child node indexes (from 1) corresponding to the columns in data.df (ignoring any grouping variable)}
\item{node.defn}{a matrix giving the parent combination}
\item{mlik}{log marginal likelihood value for each node combination. If the model cannot be fitted then NA is returned. }
\item{error.code}{if non-zero then either the root finding algorithm (glm nodes) or the maximisation algorithm (glmm nodes) terminated in an unusual way suggesting a possible unreliable result, or else the finite difference hessian estimation produced and error or warning (glmm nodes). NULL if \code{method="mle"}.}
\item{error.code.desc}{a textual description of the \code{error.code}. NULL if \code{method="mle"}}
\item{hessian.accuracy}{An estimate of the error in the final mlik value for each parent combination - this is the absolute difference between two different adaptive finite difference rules where each computes the mlik value. NULL if \code{method="mle"}}
\item{data.df}{a version of the original data (for internal use only in other functions such as \code{\link{mostprobable}}).}
\item{data.dists}{the named list of nodes distributions (for internal use only in other functions such as \code{\link{mostprobable}}).}
\item{max.parents}{the maximum number of parents (for internal use only in other functions such as \code{\link{mostprobable}}).}
\item{dag.retained}{the matrix encoding the retained arcs (for internal use only in other functions such as \code{\link{search.heuristic}}).}
\item{dag.banned}{the matrix encoding the banned arcs (for internal use only in other functions such as \code{\link{search.heuristic}}).}
\item{aic}{aic value for each node combination. If the model cannot be fitted then NaN is returned. NULL if \code{method="bayes"}.}
\item{bic}{bic value for each node combination. If the model cannot be fitted then NaN is returned. NULL if \code{method="bayes"}.}
\item{mdl}{mdl value for each node combination. If the model cannot be fitted then NaN is returned. NULL if \code{method="bayes"}.}
}

\details{
The function computes a cache of scores based on possible restrictions (maximum complexity, retained and banned arcs).

This function is very similar to \code{\link{fitAbn}} - see that help
page for details of the type of models used and in particular
\code{data.dists} specification - but rather than fit a single complete
DAG \code{buildScoreCache} iterates over all different parent
combinations for each node, creating a cache of scores.  This cache of
score could be used to select the optimal network in other function such
as \code{\link{searchHeuristic}} or \code{\link{mostprobable}}.


Two very different approaches are implemented: a Bayesian and frequentist approaches. They can be selected using the \code{method} argument.

If \code{method="bayes"}:
This function is used to calculate all individual node scores (log
marginal likelihoods).

The variable \code{which.nodes} is to allow the computation to be separated by node, for example, over different CPUs using say \code{R CMD BATCH}. This may useful and indeed likely essential with larger problems or those with random effects. Note that in this case, the results must then be combined back into a list of identical formats to that produced by an individual call to \code{buildScoreCache}, comprising of all nodes (in the same order as the columns in \code{data.df}) before sending it to any search routines. Using \code{dry.run} can be useful here.

If \code{method="mle"}:
This function is used to calculate all individual information-theoretic node scores. The possible information-theoretic based network scores
computed in \code{buildScoreCache} are the maximum likelihood (mlik, called marginal likelihood in this context as it is computed node wise),
the Akaike Information Criteria (aic), the Bayesian Information Criteria (bic) and the Minimum distance Length (mdl). The classical definitions of those metrics are given in Kratzer and Furrer (2018). This function computes a cache that can be fed into a model search algorithm.

The numerical routines used here are identical to those in \code{\link{fitAbn}} and see that help page for further details and also the quality assurance section on the \url{http://r-bayesian-networks.org} of the \pkg{abn} website for more details.

}

\references{
Lewis, F. I., and McCormick, B. J. J. (2012). "Revealing the complexity of health determinants in resource poor settings". \emph{American Journal Of Epidemiology}. doi:10.1093/aje/KWS183).

Kratzer, G., Lewis, F.I., Comin, A., Pittavino, M., and Furrer, R. (2019). "Additive Bayesian Network Modelling with the R Package abn". arXiv:1911.09006.

Kratzer, G., and Furrer, R., (2018). "Information-Theoretic Scoring Rules to Learn Additive Bayesian Network Applied to Epidemiology". arXiv:1808.01126.


Further information about \pkg{abn} can be found at: \url{http://r-bayesian-networks.org}.}

\author{Fraser Iain Lewis and Gilles Kratzer}
\seealso{\code{\link{buildScoreCache}}}
\examples{
#################################################################
## Example 1
#################################################################

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

res.c <- buildScoreCache(data.df=mydat, data.dists=mydists,
              dag.banned= ~b1|b2, dag.retained= ~g1|g2, max.parents=max.par)


## repeat but using R-INLA. The mlik's should be virtually identical.
## now build cache:
if(requireNamespace("INLA", quietly = TRUE)){  
  res.inla <- buildScoreCache(data.df=mydat, data.dists=mydists,
                dag.banned= ~b1|b2, dag.retained= ~g1|g2, max.parents=max.par,
                control=list(max.mode.error=100))

  ## comparison - very similar
  difference <- res.c$mlik - res.inla$mlik
}

## Comparison Bayes with MLE (unconstrained):
res.mle <- buildScoreCache(data.df=mydat, data.dists=mydists,
                           max.parents=3, method="mle")
res.abn <- buildScoreCache(data.df=mydat, data.dists=mydists,
                           max.parents=3, method="Bayes")
## of course different, but smame order:
plot(-res.mle$bic, res.abn$mlik)

\dontrun{
#################################################################
## Example 2 - mle with several cores
###################################################################

## Many variables, few observations
mydat <- ex0.dag.data
mydists <- as.list(rep(c("binomial", "gaussian", "poisson"), each=10))
names(mydists) <- names(mydat)

# system.time( {
# res.mle1 <- buildScoreCache(data.df=mydat, data.dists=mydists,
#                           max.parents=2, method="mle", ncores=2) })
# system.time( {
# res.mle2 <- buildScoreCache(data.df=mydat, data.dists=mydists,
#                           max.parents=2, method="mle") })


#################################################################
## Example 3 - grouped data - random effects example e.g. glmm
###################################################################

mydat <- ex3.dag.data ## this data comes with abn see ?ex3.dag.data

mydists <- list(b1="binomial", b2="binomial", b3="binomial",
            b4="binomial", b5="binomial", b6="binomial", b7="binomial",
            b8="binomial", b9="binomial", b10="binomial",b11="binomial",
            b12="binomial", b13="binomial" )
max.par <- 2

## in this example INLA is used as default since these are glmm nodes
## when running this at node-parent combination 71 the default accuracy check on the
## INLA modes is exceeded (default is a max. of 10 percent difference from
## modes estimated using internal code) and a message is given that internal code
## will be used in place of INLA's results.

# mycache <- buildScoreCache(data.df=mydat, data.dists=mydists, group.var="group",
#                         cor.vars=c("b1","b2","b3","b4","b5","b6","b7",
#                                    "b8","b9","b10","b11","b12","b13"),
#                         max.parents=max.par, which.nodes=c(1))
}
}

\keyword{models}
\concept{abn}
