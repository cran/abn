% fitabn.Rd --- 
% Author    :       Fraser Lewis
% Maintainer:       Marta Pittavino
% Last Modified By: Marta Pittavino
% Last Modified On: Gilles Kratzer
% Last modified on : 26/09/2012 
% Last modified on : 22/05/2019 GK
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

\name{fitabn}
\alias{fitabn}
\encoding{latin1}
%- Also NEED an `\alias' for EACH other topic documented here.

\title{Fit an additive Bayesian network model}

\description{Fits an additive Bayesian network to observed data and is equivalent to Bayesian or information-theoretic multi-dimensional regression modeling. Two numerical options are available in the Bayesian settings, standard Laplace approximation or else an integrated nested Laplace approximation provided via a call to the R INLA library (see www.r-inla.org - this is not hosted on CRAN).}

\usage{
      fitabn(object = NULL, dag.m = NULL, data.df = NULL, 
             data.dists = NULL, method = "bayes", group.var = NULL,
             adj.vars = NULL, cor.vars = NULL, create.graph = FALSE,
             compute.fixed = FALSE, control = list(), verbose = FALSE,
             centre = TRUE, ...)
}
  
%- maybe also `usage' for other objects documented here.
\arguments{ 
\item{object}{an object of class 'abnLearned' produced by \code{\link{mostprobable}}, \code{\link{search.heuristic}} or \code{\link{search.hillclimber}}.}
  \item{dag.m}{a matrix or a formula statement (see details) defining the network structure, a directed acyclic graph (DAG), see details for format. Note that column names and row names must be set up.}
  \item{data.df}{a data frame containing the data used for learning the network, binary variables must be declared as factors, and no missing values all allowed in any variable.}
  \item{data.dists}{a named list giving the distribution for each node in the network, see details.}
    \item{method}{should a "Bayes" or "mle" approach be used, see details.}
  \item{group.var}{only applicable for mixed models and gives the column name in \code{data.df} of the grouping variable (which must be a factor denoting group membership).}
      \item{adj.vars}{a character vector giving the column names in \code{data.df} for which the network score has to be adjusted for, see details.}
  \item{cor.vars}{a character vector giving the column names in data.df for which a mixed model should be used.}
  \item{create.graph}{create an R graph class object which enables easy plotting of dag.m using \code{plot()}, requires Rgraphviz (hosted on bioconductor rather than CRAN).}
  \item{compute.fixed}{a logical flag, set to TRUE for computation of marginal posterior distributions, see details.}
    \item{control}{a list of control parameters. See details.}
    \item{verbose}{if \code{TRUE} then provides some additional output, in particular the code used to call INLA, if applicable.}
  \item{centre}{should the observations in each Gaussian node first be standardised to mean zero and standard deviation one.}
    \item{\dots}{additional arguments passed for optimization.}
}

\details{

If \code{method="Bayes"}:

 The procedure \code{fitabn} fits an additive Bayesian network model to data where each node (variable - a column in data.df) can be either: presence/absence (Bernoulli); continuous (Gaussian); or an unbounded count (Poisson). The model comprises of a set of conditionally independent generalized linear regressions with or without random effects. Internal code is used by default for numerical estimation in nodes without random effects, and INLA is the default for nodes with random effects. This default behavior can be overridden using max.mode.error. The default is \code{max.mode.error=10}, which means that the modes estimated from INLA output must be within 10\% of those estimated using internal code. Otherwise, the internal code is used rather than INLA. To force the use of INLA on all nodes, use max.mode.error=100, which then ignores this check, to force the use of internal code then use \code{max.mode.error=0}. For the numerical reliability and perform of abn see \url{http://r-bayesian-networks.org}. Generally speaking, INLA can be swift and accurate, but in several cases, it can perform very poorly and so some care is required (which is why there is an internal check on the modes). Binary variables must be declared as factors with two levels, and the argument \code{data.dists} must be a list with named arguments, one for each of the variables in \code{data.df} (except a grouping variable - if present), where each entry is either "poisson","binomial", or "gaussian", see examples below. The "poisson" and "binomial" distributions use log and logit link functions, respectively. Note that "binomial" here actually means only binary, one Bernoulli trial per row in \code{data.df}. 

If the data are grouped into correlated blocks - wherein a standard regression context a mixed model might be used - then a network comprising of one or more nodes where a generalized linear mixed model is used (but limited to only a single random effect). This is achieved by specifying parameters \code{group.var} and \code{cor.vars}. Where the former defines the group membership variable, which should be a factor indicating which observations belong to the same grouping. The parameter \code{cor.vars} is a character vector that contains the names of the nodes for which a mixed model should be used. For example, in some problems, it may be appropriate for all variables (except \code{group.var}) in data.df to be parametrized as a mixed model while in others it may only be a single variable for which grouping adjustment is required (as the remainder of variables are covariates measured at group level).   

In the network structure definition, \code{dag.m}, each row represents a node in the network, and the columns in each row define the parents for that particular node, see the example below for the specific format. The \code{dag.m} can be provided using a formula statement (similar to GLM). A typical formula is \code{ ~ node1|parent1:parent2 + node2:node3|parent3}. The formula statement have to start with \code{~}. In this example, node1 has two parents (parent1 and parent2). node2 and node3 have the same parent3. The parents names have to exactly match those given in \code{data.df}. \code{:} is the separator between either children or parents, \code{|} separates children (left side) and parents (right side), \code{+} separates terms, \code{.} replaces all the variables in \code{data.df}.

If \code{compute.fixed=TRUE} then the marginal posterior distributions for all parameters are computed. Note the current algorithm used to determine the evaluation grid is rather crude and may need to be manually refined using \code{variate.vec} (one parameter at a time) for publication-quality density estimates. Note that a manual grid can only be used with internal code and not INLA (which uses its own grid). The end points are defined as where the value of the marginal density drops below a given threshold \code{pdf.min}. 

If \code{create.graph=TRUE} then the model definition matrix in \code{dag.m} is used to create an R graph object (of type graphAM-class). See ?"graph-class" for details and the Rgraphviz documentation (which is extensive). The main purpose of this is to allow easy visualization of the DAG structure via the \code{graphviz} library. A graph plot can easily be created by calling the method \code{plot} on this object (see example below). Note, however, that the fonts and choice of scaling used here may be far less visually optimal than using \code{graphviz} direct (e.g., via \code{\link{tographviz}}) for publication-quality output. Also, re-scaling the plotting window may not result in a callback to re-optimize the visual position of the nodes and edges, and so if the window is re-sized, then re-run the plot command to re-draw to the new scale.   

When estimating the log marginal likelihood in models with random effects (using internal code rather than INLA), an attempt is made to minimize the error by comparing the estimates given between a 3pt and 5pt rule when estimating the Hessian in the Laplace approximation. The modes used in each case are identical. The first derivatives are computed using gsl's adaptive finite difference function, and this is embedding inside the standard 3pt and 5pt rules for the second derivatives. In all cases, a central difference approximation is tried first with a forward difference being a fall back (as the precision parameters are strictly positive). The error is minimized through choosing an optimal step size using gsl's Nelder-Mead optimization, and if this fails, (e.g., is larger than \code{max.hessian.error}) then the Brent-Dekker root bracketing method is used as a fallback. If the error cannot be reduced to below \code{max.hessian.error}, then the step size, which gave the lowest error during the searches (across potentially many different initial bracket choices), is used for the final Hessian evaluations in the Laplace approximation.


If \code{method="mle"}:

 The procedure \code{fitabn} with the argument \code{method= "mle"} fits an additive Bayesian network model to data where each node (variable - a column in data.df) can be either: presence/absence (Bernoulli); continuous (Gaussian); an unbounded count (Poisson); or a discrete variable (Multinomial). The model comprises of a set of conditionally independent generalized linear regressions with or without adjustment. 
 
Binary and discrete variables must be declared as factors and the argument \code{data.dists} must be a list with named arguments, one for each of the variables in \code{data.df}, where each entry is either "poisson","binomial", "multinomial" or "gaussian", see examples below. The "poisson" and "binomial" distributions use log and logit link functions, respectively. Note that "binomial" here actually means only binary, one Bernoulli trial per row in data.df.  

In the context of \code{fitabn} adjustment means that irrespective to the adjacency matrix the adjustment variable set (\code{adj.vars}) will be add as covariate to every node defined by \code{cor.vars}. If \code{cor.vars} is NULL then adjustment is over all variables in the \code{data.df}.

In the network structure definition, \code{dag.m}, each row represents a node in the network, and the columns in each row define the parents for that particular node, see the example below for the specific format. The \code{dag.m} can be provided using a formula statement (similar to GLM). A typical formula is \code{ ~ node1|parent1:parent2 + node2:node3|parent3}. The formula statement have to start with \code{~}. In this example, node1 has two parents (parent1 and parent2). node2 and node3 have the same parent3. The parents names have to exactly match those given in \code{data.df}. \code{:} is the separator between either children or parents, \code{|} separates children (left side) and parents (right side), \code{+} separates terms, \code{.} replaces all the variables in \code{data.df}.

The Information-theoretic based network scores used in \code{fitabn} with argument \code{method="mle"} are the maximum likelihood (mlik, called marginal likelihood in this context as it is computed node wise), the Akaike Information Criteria (aic), the Bayesian Information Criteria (bic) and the Minimum distance Length (mdl). The classical definitions of those metrics are given in Kratzer and Furrer (2018).

The numerical routine is based on an iterative scheme to estimate the regression coefficients. The Iterative Reweighed Least Square (IRLS) programmed using Rcpp/RcppArmadrillo. One hard coded feature of \code{fitabn} with argument \code{method="mle"} is a conditional use of a bias reduced binomial regression when a classical Generalized Linear Model (GLM) fails to estimate the maximum likelihood of the given model accurately. Additionally, a QR decomposition is performed to check for rank deficiency. If the model is rank deficient and the BR GLM fails to estimate it, then predictors are sequentially removed. This feature aims at better estimating network scores when data sparsity is present.

A special care should be taken when interpreting or even displaying p-values computed with \code{fitabn}. Indeed, the full model is already selected using goodness of fit metrics based on the (same) full dataset. 

The control argument is a list that can supply any of the following components: 

\describe{ 
  \item{mean}{the prior mean for all the Gaussian additive terms for each node.}
  \item{prec}{the prior precision for all the Gaussian additive terms for each node.}
  \item{loggam.shape}{the shape parameter in the Gamma distributed prior for the precision in any Gaussian nodes, also used for group-level precision is applicable.}
  \item{loggam.inv.scale}{the inverse scale parameter in the Gamma distributed prior for the precision in any Gaussian nodes, also used for group-level precision, is applicable. }
  \item{max.mode.error}{if the estimated modes from INLA differ by a factor of max.mode.error or more from those computed internally, then results from INLA are replaced by those computed internally. To force INLA always to be used, then max.mode.error=100, to force INLA never to be used max.mod.error=0. See details.}  
  \item{max.iters}{total number of iterations allowed when estimating the modes in Laplace approximation}
  \item{epsabs}{absolute error when estimating the modes in Laplace approximation for models with no random effects.}
  \item{error.verbose}{logical, additional output in the case of errors occurring in the optimization}
  \item{epsabs.inner}{absolute error in the maximization step in the (nested) Laplace approximation for each random effect term}
  \item{max.iters.inner}{total number of iterations in the maximization step in the nested Laplace approximation}
  \item{finite.step.size}{suggested step length used in finite difference estimation of the derivatives for the (outer) Laplace approximation when estimating modes}
  \item{hessian.params}{a numeric vector giving parameters for the adaptive algorithm, which determines the optimal step size in the finite-difference estimation of the Hessian. First entry is the initial guess, second entry absolute error}
\item{max.iters.hessian}{integer, maximum number of iterations to use when determining an optimal finite difference approximation (Nelder-Mead)}
\item{max.hessian.error}{if the estimated log marginal likelihood when using an adaptive 5pt finite-difference rule for the Hessian differs by more than max.hessian.error from when using an adaptive 3pt rule then continue to minimize the local error by switching to the Brent-Dekker root bracketing method, see details} 
\item{factor.brent}{if using Brent-Dekker root bracketing method then define the outer most interval end points as the best estimate of h (stepsize) from the Nelder-Mead as (h/factor.brent,h*factor.brent)}
\item{maxiters.hessian.brent}{maximum number of iterations allowed in the Brent-Dekker method}
\item{num.intervals.brent}{the number of initial different bracket segments to try in the Brent-Dekker method}
\item{min.pdf}{the value of the posterior density function below which we stop the estimation only used when computing marginals, see details.} 
\item{n.grid}{recompute density on an equally spaced grid with \code{n.grid} points.}
\item{std.area}{logical, should the area under the estimated posterior density be standardized to exactly one, useful for error checking.} 
\item{marginal.quantiles}{vector giving quantiles at which to compute the posterior marginal distribution at.}
\item{max.grid.iter}{gives number of grid points to estimate posterior density at when not explicitly specifying a grid used to avoid excessively long computation.}
\item{marginal.node}{used in conjunction with \code{marginal.param} to allow bespoke estimate of a marginal density over a specific grid. value from 1 to the number of nodes.} 
 \item{marginal.param}{used in conjunction with \code{marginal.node}. value of 1 is for intercept, see modes entry in results for the appropriate number.}
\item{variate.vec}{a vector containing the places to evaluate the posterior marginal density, must be supplied if \code{marginal.node} is not null}
    \item{maxit}{integer given the maximum number of run for estimating network scores using an Iterative Reweighed Least Square algorithm.}
  \item{tol}{real number giving the minimal tolerance expected to terminate the Iterative Reweighed Least Square algorithm to estimate network score.}
  \item{seed}{a non-negative integer which sets the seed.}
}

}

\value{An object of class \code{abnFit}. A named list. One entry for each of the variables in \code{data.df} (excluding the grouping variable, if present) which contains an estimate of the log marginal likelihood for that individual node. An entry "mlik" which is the total log marginal likelihood for the full ABN model. A vector of \code{error.codes} - non-zero if a numerical error or warning occurred, and a vector error.code.desc giving a text description of the error. A list \code{modes}, which contains all the mode estimates for each parameter at each node. A vector called Hessian accuracy, which is the estimated local error in the log marginal likelihood for each node.  If \code{compute.fixed=TRUE} then a list entry called \code{marginals} which contains a named entry for every parameter in the ABN and each entry in this list is a two-column matrix where the first column is the value of the marginal parameter, say x, and the second column is the respective density value, pdf(x). Also, a list called \code{marginal.quantiles} is produced, giving the quantiles for each marginal posterior distribution. If \code{create.graph=TRUE} then an additional entry \code{graph} which is of type class \code{graphAM}-class (\pkg{Rgraphviz} is created. 

}

\references{Lewis, F. I., and McCormick, B. J. J. (2012). Revealing the complexity of health determinants in resource poor settings. \emph{American Journal Of Epidemiology}. DOI:10.1093/aje/KWS183.  

Kratzer, G., Lewis, F.I., Comin, A., Pittavino, M. and Furrer, R. (2019). "Additive Bayesian Network Modelling with the R Package abn". arXiv preprint arXiv:1911.09006.

Kratzer, G., and Furrer, R., 2018. Information-Theoretic Scoring Rules to Learn Additive Bayesian Network Applied to Epidemiology. Preprint; Arxiv: stat.ML/1808.01126. 


Further information about \pkg{abn} can be found at:\cr
\url{http://r-bayesian-networks.org}}

\author{Fraser Iain Lewis and Gilles Kratzer}

\seealso{\code{\link{buildscorecache}}}

\examples{
\dontrun{
## Built-in dataset with a subset of cols  
mydat <- ex0.dag.data[,c("b1","b2","b3","g1","b4","p2","p4")] 


## setup distribution list for each node
mydists <- list(b1="binomial", b2="binomial", b3="binomial", g1="gaussian",
              b4="binomial", p2="poisson", p4="poisson")

## Null model - all independent variables
mydag.empty <- matrix(0, nrow=7, ncol=7)
colnames(mydag.empty) <- rownames(mydag.empty) <- names(mydat)

## Now fit the model to calculate its goodness-of-fit
myres <- fitabn(dag.m = mydag.empty, data.df = mydat, data.dists = mydists)

## Log-marginal likelihood goodness-of-fit for complete DAG
print(myres$mlik) 

## Now repeat but include some dependencies first
mydag <- mydag.empty 
mydag["b1","b2"] <- 1 # b1<-b2 
mydag["b2","p4"] <- 1 # b2<-p4
mydag["b2","g1"] <- 1 # b2<-g1
mydag["g1","p2"] <- 1 # g1<-p2
mydag["b3","g1"] <- 1 # b3<-g1
mydag["b4","b1"] <- 1 # b4<-b1
mydag["p4","g1"] <- 1 # p4<-g1

# fit  using the formula statement 
#  including the creation of the graph of the DAG via Rgraphviz 
myres <- fitabn(dag.m = ~b1|b2+b2|p4+g1+g1|p2+b3|g1+b4|b1+p4|g1,
                data.df = mydat, data.dists = mydists,
                create.graph = TRUE)

## Or equivalentelly using the formula statement
# myres <- fitabn(dag.m = mydag, data.df = mydat, data.dists = mydists,
#              create.graph = TRUE)

print(myres$mlik) ## a much weaker fit than full independence DAG

plotabn(dag.m = mydag, data.dists = mydists, fitted.values.abn = myres$modes)

## A simple plot of some posterior densities the algorithm which chooses 
## density points is very simple any may be rather sparse so also recompute 
## the density over an equally spaced grid of 50 points between the two 
## end points which had at f=min.pdf
## max.mode.error = 0 foces to use the internal c code
myres.c <- fitabn(dag.m = mydag, data.df = mydat, data.dists = mydists, 
                  compute.fixed = TRUE, n.grid = 50, max.mode.error = 0)

print(names(myres.c$marginals)) ## gives all the different parameter names

## Repeat but use INLA for the numerics using max.mode.error=100
## as using internal code is the default here rather than INLA 
myres.inla <- fitabn(dag.m = mydag, data.df = mydat, data.dists = mydists,
                     compute.fixed = TRUE, n.grid = 50, max.mode.error = 100)

## Plot posterior densities
par(mfrow=c(2,2), mai=c(.7,.7,.2,.1))
plot(myres.c$marginals$b1[["b1|(Intercept)"]], type="l", xlab="b1|(Intercept)")
lines(myres.inla$marginals$b1[["b1|(Intercept)"]], col="blue")
plot(myres.c$marginals$b2[["b2|p4"]], type="l", xlab="b2|p4")
lines(myres.inla$marginals$b2[["b2|p4"]], col="blue")
plot(myres.c$marginals$g1[["g1|precision"]], type="l", xlab="g1|precision")
lines(myres.inla$marginals$g1[["g1|precision"]], col="blue")
plot(myres.c$marginals$b4[["b4|b1"]], type="l", xlab="b4|b1")
lines(myres.inla$marginals$b4[["b4|b1"]], col="blue")

## An elementary mixed model example using built-in data specify DAG, 
## only two variables using a subset of variables from ex3.dag.data
## both variables are assumed to need (separate) adjustment for the 
## group variable, i.e., a binomial GLMM at each node


mydists <- list(b1="binomial",  b2="binomial")

## Compute marginal likelihood - use internal code via max.mode.error=0
## as using INLA is the default here.
## Model where b1 <- b2
myres.c <- fitabn(dag.m = ~b1|b2, data.df = ex3.dag.data[,c(1,2,14)], data.dists = mydists,
                  group.var = "group", cor.vars = c("b1","b2"),
                  max.mode.error=0)
print(myres.c) ## show all the output 

## compare mode for node b1 with glmer(), lme4::glmer is automatically attached.
m1 <- glmer(b1 ~ 1 + b2 + (1|group), 
            family = "binomial", data = ex3.dag.data[,c(1,2,14)])

print(myres.c$modes$b1) ## almost identical to lme4 n.b. glmer() gives variance
##                         fitabn gives precision=1/variance

## Compare with INLA estimate
myres.inla <- fitabn(dag.m = ~b1|b2,data.df=ex3.dag.data[,c(1,2,14)],
                     data.dists = mydists, group.var = "group", cor.vars = c("b1","b2"),
                     compute.fixed = FALSE, max.mode.error = 100)

## Compare log marginal likelihoods for each node and total DAG:
cbind(INLA = unlist(myres.inla[1:3]), C = unlist(myres.c[1:3]), 
             Delta = unlist(myres.inla[1:3]) - unlist(myres.c[1:3]))

## Now for marginals - INLA is strongly preferable for estimating marginals for nodes 
## with random effects as it is far faster, but may not be reliable
## see http://r-bayesian-networks.org

## INLA's estimates of the marginals, using high n.grid=500
## as this makes the plots smoother - see below.
## myres.inla <- fitabn(dag.m = ~b1|b2, data.df = ex3.dag.data[,c(1,2,14)], 
#                      data.dists = mydists,
#                      group.var = "group", cor.vars = c("b1", "b2"),
#                      compute.fixed = TRUE, n.grid = 500, 
#                      control = list(max.mode.error = 100, max.hessian.error = 10E-02))

## this is NOT recommended - marginal density estimation using fitabn in mixed models
## is really just for diagnostic purposes, better to use fitabn.inla() here
## but here goes...be patient
# myres.c <- fitabn(dag.m = ~b1|b2, data.df = ex3.dag.data[,c(1,2,14)], data.dists = mydists,
#                   group.var = "group", cor.vars=c("b1", "b2"), compute.fixed = TRUE,
#                   control = list(max.mode.error = 0, max.hessian.error = 10E-02))

## compare marginals between internal and INLA.   
# par(mfrow=c(2,3))
## 5 parameters - two intercepts, one slope, two group level precisions
# plot(myres.inla$marginals$b1[[1]], type="l", col="blue")
# lines(myres.c$marginals$b1[[1]], col="brown", lwd=2)
# plot(myres.inla$marginals$b1[[2]], type="l", col="blue")
# lines(myres.c$marginals$b1[[2]], col="brown", lwd=2)
## the precision of group-level random effects
# plot(myres.inla$marginals$b1[[3]],type="l", col="blue", xlim=c(0,2))
# lines(myres.c$marginals$b1[[3]],col="brown",lwd=2)
# plot(myres.inla$marginals$b2[[1]],type="l", col="blue")
# lines(myres.c$marginals$b2[[1]],col="brown",lwd=2)
# plot(myres.inla$marginals$b2[[1]], type="l", col="blue")
# lines(myres.c$marginals$b2[[1]], col="brown", lwd=2)
## the precision of group-level random effects
# plot(myres.inla$marginals$b2[[2]], type="l", col="blue", xlim=c(0,2))
# lines(myres.c$marginals$b2[[2]], col="brown", lwd=2)

### these are very similar although not exactly identical

## use internal code but only to compute a single parameter over a specified grid
## This can be necessary if the simple auto grid finding functions does a poor job

#myres.c <- fitabn(dag.m = ~b1|b2, data.df = ex3.dag.data[,c(1,2,14)], data.dists = mydists,
#                  group.var = "group", cor.vars = c("b1", "b2"), 
#                  centre = FALSE, compute.fixed = TRUE,
#                  marginal.node = 1, marginal.param = 3,## precision term in node 1
#                  variate.vec = seq(0.05, 1.5, len=25), max.hessian.error = 10E-02)

#par(mfrow=c(1,2))
#plot(myres.c$marginals[[1]], type="l", col="blue")## still fairly sparse
## An easy way is to use spline to fill in the density without recomputing other
## points provided the original grid is not too sparse.
#plot(spline(myres.c$marginals[[1]], n=100), type="b", col="brown")

## -----------------------------------------------------------------------------------
## This function contains an MLE implementation accessible through a method parameter
## use built-in simulated data set
## -----------------------------------------------------------------------------------

mydat <- ex0.dag.data[,c("b1","b2","b3","g1","b4","p2","p4")] ## take a subset of cols

## setup distribution list for each node
mydists <- list(b1="binomial", b2="binomial", b3="binomial",
                g1="gaussian", b4="binomial", p2="poisson", p4="poisson")


## now fit the model to calculate its goodness of fit
myres.mle <- fitabn(dag.m = ~b1|b2+b2|p4+g1+g1|p2+b3|g1+b4|b1+p4|g1, 
                    data.df = mydat, data.dists = mydists, method = "mle")

myres.bayes <- fitabn(dag.m = ~b1|b2+b2|p4+g1+g1|p2+b3|g1+b4|b1+p4|g1, 
                      data.df = mydat, data.dists = mydists, method = "bayes")

## print the output
## MLE
print(myres.mle)

#Bayes
print(myres.bayes)

## plot the model with parameter estimates
plotabn(dag.m = mydag, data.dists = mydists, fitted.values.abn.mle = myres.bayes$modes)
}

}

\keyword{models}
\concept{abn}