% buildscorecache.Rd --- 
% Author           : Fraser Lewis & Gilles Kratzer
% Created on :       03/08/2012 
% Last modification : 01/12/2016
% Last modification : 01/03/2017
% Last modification : 21/05/2019 (GK)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

\name{buildscorecache}

\alias{buildscorecache}

\title{Build a cache of goodness of fit metrics for each node in a DAG, possibly subject to user-defined restrictions}

\description{
  Iterates over all valid parent combinations - subject to ban, retain, and \code{max.parent} limits - for each node, or a subset of nodes, and computes a cache of log marginal likelihoods. This cache can then be used in different DAG structural search algorithms. 
}

\usage{
buildscorecache(data.df = NULL, data.dists = NULL, method = "bayes",
   group.var = NULL, adj.vars = NULL, cor.vars = NULL, dag.banned = NULL,
   dag.retained = NULL, max.parents = NULL, which.nodes=NULL,
   defn.res = NULL, centre = TRUE, dry.run = FALSE,
   control = list(), verbose = FALSE, ...)
}
\arguments{
  \item{data.df}{a data frame containing the data used for learning each node, binary variables must be declared as factors.}
  \item{data.dists}{a named list giving the distribution for each node in the network, see \sQuote{Details}.} 
  \item{method}{should a "Bayes" or "mle" approach be used, see \sQuote{Details}.}
  \item{group.var}{only applicable for nodes to be fitted as a mixed model and gives the column name in \code{data.df} of the grouping variable which must be a factor denoting group membership.}
    \item{adj.vars}{a character vector giving the column names in \code{data.df} for which the network score has to be adjusted for, see \sQuote{Details}.}
  \item{cor.vars}{a character vector giving the column names in \code{data.df} for which a mixed model should be used to adjust for within group correlation or pure adjustment.}
  \item{dag.banned}{a matrix or a formula statement (see \sQuote{Details} for format) defining which arcs are not permitted - banned - see \sQuote{Details} for format. Note that colnames and rownames must be set, otherwise same row/column names as data.df will be assumed. If set as NULL an empty matrix is assumed.}
  \item{dag.retained}{a matrix or a formula statement (see \sQuote{Details} for format) defining which arcs are must be retained in any model search, see \sQuote{Details} for format. Note that colnames and rownames must be set, otherwise same row/column names as data.df will be assumed. If set as NULL an empty matrix is assumed.}
  \item{max.parents}{a constant or named list giving the maximum number of parents allowed, the list version allows this to vary per node.}
  \item{which.nodes}{a vector giving the column indices of the variables to be included, if ignored all variables are included.}
  \item{defn.res}{an optional user-supplied list of child and parent combinations, see \sQuote{Details}.}
    \item{centre}{should the observations in each Gaussian node first be standardized to mean zero and standard deviation one, defaults to TRUE.}
  \item{dry.run}{if TRUE then a list of the child nodes and parent combinations are returned but without estimation of node scores (log marginal likelihoods).}
  \item{control}{a list of control parameters. See \sQuote{Details}.}
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
\item{aic}{aic value for each node combination. If the model cannot be fitted then NaN is returned. NULL if \code{method="Bayes"}.}
\item{bic}{bic value for each node combination. If the model cannot be fitted then NaN is returned. NULL if \code{method="Bayes"}.}
\item{mdl}{mdl value for each node combination. If the model cannot be fitted then NaN is returned. NULL if \code{method="Bayes"}.}
}

\details{
The function computes a cache of scores based on possible restrictions (maximum complexity, retained and banned arcs). This cache of score could be used to select the optimal network in other function such as \code{\link{searchHeuristic}} or \code{\link{mostprobable}}. Two very different approaches are implemented: a Bayesian and frequentist approaches. They can be selected using the \code{method} argument. 

If \code{method="Bayes"}:

This function is used to calculate all individual node scores (log marginal likelihoods). This cache can then be fed into a model search algorithm. This function is very similar to \code{\link{fitabn}} - see that help page for details of the type of models used and in particular \code{data.dists} specification - but rather than fit a single complete DAG \code{buildscorecache} iterates over all different parent combinations for each node.  There are three ways to customise the parent combinations through giving a matrix which contains arcs which are not allowed (banned), a matrix which contains arcs which must always be included (retained) and also a general complexity limit which restricts the maximum number of arcs allowed to terminate at a node (its number of parents), where this can differ from node to node. In these matrices, \code{dag.banned} and \code{dag.retained}, each row represents a node in the network, and the columns in each row define the parents for that particular node, see the example below for the specific format. If these are not supplied, they are assumed to be empty matrices, i.e., no arcs banned or retained. Note that only rudimentary consistency checking is done here, and some care should be taken not to provide conflicting restrictions in the ban and retain matrices.

The variable \code{which.nodes} is to allow the computation to be separated by node, for example, over different CPUs using say \code{R CMD BATCH}. This may useful and indeed likely essential with larger problems or those with random effects. Note that in this case, the results must then be combined back into a list of identical formats to that produced by an individual call to \code{buildscorecache}, comprising of all nodes (in the same order as the columns in \code{data.df}) before sending it to any search routines. Using \code{dry.run} can be useful here.  

If \code{method="mle"}:

This function is used to calculate all individual Information-Theoretic node scores. The possible Information-theoretic based network scores computed in \code{buildscorecache} are the maximum likelihood (mlik, called marginal likelihood in this context as it is computed node wise), the Akaike Information Criteria (aic), the Bayesian Information Criteria (bic) and the Minimum distance Length (mdl). The classical definitions of those metrics are given in Kratzer and Furrer (2018). This function computes a cache that can be fed into a model search algorithm. This function is very similar to \code{\link{fitabn}} with \code{method="mle"} - see that help page for details of the type of models used and in particular \code{data.dists} specification - but rather than fit a single complete DAG \code{buildscorecache} iterates over all admissible parent combinations for each node.  There are three ways to customise the parent combinations through giving a matrix which contains arcs which are not allowed (banned), a matrix which contains arcs which must always be included (retained) and also a general complexity limit which restricts the maximum number of arcs allowed to terminate at a node (its number of parents). In these matrices, \code{dag.banned} and \code{dag.retained}, each row represents a node in the network, and the columns in each row define the parents for that particular node, see the example below for the specific format. If these are not supplied they are assumed to be empty matrices, i.e., no arcs banned or retained. Note that only rudimentary consistency checking is done here, and some care should be taken not to provide conflicting restrictions in the ban and retain matrices.

The control argument is a list that can supply any of the following components: 

\describe{ 
 \item{max.mode.error}{if the estimated modes from INLA differ by a factor of max.mode.error or more from those computed internally, then results from INLA are replaced by those computed internally. To force INLA always to be used, then max.mode.error=100, to force INLA never to be used max.mod.error=0. See \sQuote{Details}.}
  \item{mean}{the prior mean for all the Gaussian additive terms for each node}
  \item{prec}{the prior precision for all the Gaussian additive term for each node}
  \item{loggam.shape}{the shape parameter in the Gamma distribution prior for the precision in a Gaussian node}
  \item{loggam.inv.scale}{the inverse scale parameter in the Gamma distribution prior for the precision in a Gaussian node}
   \item{max.iters}{total number of iterations allowed when estimating the modes in Laplace approximation}
  \item{epsabs}{absolute error when estimating the modes in Laplace approximation for models with no random effects.}
  \item{error.verbose}{logical, additional output in the case of errors occurring in the optimization}
   \item{epsabs.inner}{absolute error in the maximization step in the (nested) Laplace approximation for each random effect term}
   \item{max.iters.inner}{total number of iterations in the maximization step in the nested Laplace approximation}
  \item{finite.step.size}{suggested step length used in finite difference estimation of the derivatives for the (outer) Laplace approximation when estimating modes}
  \item{hessian.params}{a numeric vector giving parameters for the adaptive algorithm, which determines the optimal stepsize in the finite-difference estimation of the hessian. First entry is the initial guess, second entry absolute error}
\item{max.iters.hessian}{integer, maximum number of iterations to use when determining an optimal finite difference approximation (Nelder-Mead)}
\item{max.hessian.error}{if the estimated log marginal likelihood when using an adaptive 5pt finite-difference rule for the Hessian differs by more than max.hessian.error from when using an adaptive 3pt rule then continue to minimize the local error by switching to the Brent-Dekker root bracketing method, see \sQuote{Details}} 
\item{factor.brent}{if using Brent-Dekker root bracketing method then define the outer most interval end points as the best estimate of h (stepsize) from the Nelder-Mead as (h/factor.brent,h*factor.brent)}
\item{maxiters.hessian.brent}{maximum number of iterations allowed in the Brent-Dekker method}
\item{num.intervals.brent}{the number of initial different bracket segments to try in the Brent-Dekker method}
  \item{maxit}{integer given the maximum number of run for estimating network scores using an Iterative Reweighed Least Square algorithm.}
  \item{tol}{real number giving the minimal tolerance expected to terminate the Iterative Reweighed Least Square algorithm to estimate network score.}
  \item{seed}{a non-negative integer which sets the seed.}
}


The numerical routines used here are identical to those in \code{\link{fitabn}} and see that help page for further details and also the quality assurance section on the \url{http://r-bayesian-networks.org} of the abn website for more details.
 
}

\references{
Lewis, F. I., and McCormick, B. J. J. (2012). "Revealing the complexity of health determinants in resource poor settings". \emph{American Journal Of Epidemiology}. DOI:10.1093/aje/KWS183).  

Kratzer, G., Lewis, F.I., Comin, A., Pittavino, M. and Furrer, R. (2019). "Additive Bayesian Network Modelling with the R Package abn". arXiv preprint arXiv:1911.09006.

Kratzer, G., and Furrer, R., (2018). "Information-Theoretic Scoring Rules to Learn Additive Bayesian Network Applied to Epidemiology". Preprint; Arxiv: stat.ML/1808.01126. 


  Further information about \pkg{abn} can be found at:\cr
  \url{http://r-bayesian-networks.org}}

\author{Fraser Iain Lewis and Gilles Kratzer}

\examples{
\dontrun{
#################################################################
## Example 1
#################################################################

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

res.c <- buildscorecache(data.df=mydat, data.dists=mydists,
              dag.banned= ~b1|b2, dag.retained= ~g1|g2, max.parents=max.par)


## repeat but using R-INLA. The mlik's should be virtually identical.
## now build cache
res.inla <- buildscorecache(data.df=mydat, data.dists=mydists,
              dag.banned= ~b1|b2, dag.retained= ~g1|g2, max.parents=max.par,
              max.mode.error=100)

## plot comparison - very similar
plot(res.c$mlik, res.inla$mlik, pch="+") 
abline(0, 1)


#################################################################
## Example 2 - much bigger problem using glm - may take a while
#################################################################

mydat <- ex2.dag.data ## this data comes with abn see ?ex2.dag.data

## setup distribution list for each node
mydists <- list(b1="binomial", g1="gaussian", p1="poisson",
                b2="binomial", g2="gaussian", p2="poisson",
                b3="binomial", g3="gaussian", p3="poisson",
                b4="binomial", g4="gaussian", p4="poisson",
                b5="binomial", g5="gaussian", p5="poisson",
                b6="binomial", g6="gaussian", p6="poisson")


## no explicit ban or retain restrictions set so do not need to supply ban
##  or retain matrices

## now build cache using internal code just for nodes 1,2 and 3
## e.g. "b1", "p1" and "g1" 
mycache.c <- buildscorecache(data.df=mydat, data.dists=mydists,
                         max.parents=2, which.nodes=c(1:3))

###################################################################
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

# mycache <- buildscorecache(data.df=mydat, data.dists=mydists, group.var="group",
#                         cor.vars=c("b1","b2","b3","b4","b5","b6","b7",
#                                    "b8","b9","b10","b11","b12","b13"),
#                         max.parents=max.par, which.nodes=c(1))
}
mydat <- ex0.dag.data[,c("b1","b2","g1","g2","b3","g3")] ## take a subset of cols

## setup distribution list for each node
mydists <- list(b1="binomial", b2="binomial", g1="gaussian",
              g2="gaussian", b3="binomial", g3="gaussian")
             
## now build cache of scores (goodness of fits for each node)
res.mle <- buildscorecache(data.df=mydat, data.dists=mydists,
                           max.parents=3, method="mle")
res.abn <- buildscorecache(data.df=mydat, data.dists=mydists,
                           max.parents=3, method="Bayes")

#plot(-res.mle$bic, res.abn$mlik)
}

\keyword{models}
\concept{abn}