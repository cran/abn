% fitabn_mle.Rd --- 
% Author    :       Gilles Kratzer
% Last Modified By: 
% Last Modified On: 
% Last modified on : 18/07/2018 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

\name{fitabn.mle}
\alias{fitabn.mle}
\encoding{latin1}
%- Also NEED an `\alias' for EACH other topic documented here.

\title{Fit an additive Bayesian network model based on maximum likelihood estimation.}

\description{Fits an additive Bayesian network to observed data and is equivalent to classical generalized Linear multi-dimensional regression modelling. The estimation is done through an Iterative Reweighed Least Square scheme.}

\usage{

fitabn.mle(dag.m = NULL,
       data.df = NULL, 
       data.dists = NULL,
       adj.vars = NULL,
       cor.vars = NULL,
       centre = TRUE,
       maxit = 100, 
       tol = 10^-11)
}
  
%- maybe also `usage' for other objects documented here.
\arguments{ 
  \item{dag.m}{a matrix or a formula statement defining the network structure, a directed acyclic graph (DAG), see details for format. Note that colnames and rownames must be set.}
  \item{data.df}{a data frame containing the data used for learning the network, binary variables must be declared as factors and no missing values all allowed in any variable.}
  \item{data.dists}{a named list giving the distribution for each node in the network, see details.}
    \item{adj.vars}{a character vector giving the column names in \code{data.df} for which the network score has to be adjusted for, see details.}
  \item{cor.vars}{a character vector giving the column names in \code{data.df} for which adjustment should be used.}
  \item{centre}{logical variable, should the observations in each Gaussian node first be standardised to mean zero and standard deviation one, defaults is TRUE.}
    \item{maxit}{integer given the maximum number of run for estimating network scores using an Iterative Reweighed Least Square algorithm.}
  \item{tol}{real number giving the minimal tolerance expected to terminate the Iterative Reweighed Least Square algorithm to estimate network score.}
}

\details{
 The procedure \code{fitabn.mle} fits an additive Bayesian network model to data where each node (variable - a column in data.df) can be either: presence/absence (Bernoulli); continuous (Gaussian); an unbounded count (Poisson); or a discrete variable (Multinomial). The model comprises of a set of conditionally independent generalized linear regressions with or without adjustment. 
 
Binary and discrete variables must be declared as factors, and the argument data.dists must be a list with named arguments, one for each of the variables in \code{data.df}, where each entry is either "poisson","binomial", "multinomial" or "gaussian", see examples below. The "poisson" and "binomial" distributions use log and logit link functions respectively. Note that "binomial" here actually means only binary, one Bernoulli trial per row in data.df.  

In the context of \code{fitabn.mle} adjustment means that irrespective to the adjacency matrix the adjustment variable set (\code{adj.vars}) will be add as covariate to every node defined by \code{cor.vars}. If \code{cor.vars} is NULL then adjustment is over all variables in the \code{data.df}.

In the network structure definition, \code{dag.m}, each row represents a node in the network, and the columns in each row define the parents for that particular node, see the example below for the specific format. The \code{dag.m} can be provided using a formula statement (similar to glm). A typical formula is \code{ ~ node1|parent1:parent2 + node2:node3|parent3}. The formula statement have to start with \code{~}. In this example, node1 has two parents (parent1 and parent2). node2 and node3 have the same parent3. The parents names have to exactly match those given in \code{data.df}. \code{:} is the separator between either children or parents, \code{|} separates children (left side) and parents (right side), \code{+} separates terms, \code{.} replaces all the variables in \code{data.df}.

The Information-theoretic based network scores used in \code{fitabn.mle} are the maximum likelihood (mlik, called marginal likelihood in this context as it is computed node wise), the Akaike Information Criteria (aic), the Bayesian Information Criteria (bic) and the Minimum distance Length (mdl). The classical definitions of those metrics are given in Kratzer and Furrer (2018).

The numerical routine are based on an iterative scheme to estimate the regression coefficients. The Iterative Reweighed Least Square (IRLS) programmed using Rcpp/RcppArmadrillo. One hard coded feature of \code{fitabn.mle} is a conditional use of a bias reduced binomial regression when a classical glm fails to accurately estimates the maximum likelihood of the given model. Additionally a QR decomposition is performed to check for rank deficiency. If the model is rank deficient and the BR glm fails to estimate it, then predictors are sequentially removed. This feature aims at better estimating network scores when data sparsity is present.

A special care should be taken when intrepreting or even displaing p-values computed with \code{fitabn.mle}. Indeed, the full model is already selected using goodness of fit metrics based on the (same) full dataset.
}

\value{A named list. One entry for each of the variables in \code{data.df} (excluding the adjusting variable, if present) which contains an estimate of the log marginal likelihood for that individual node. An entry "mlik" which is the total log marginal likelihood for the full ABN model. Same entries for aic, bic and mdl. The regression coefficients, standard error and P-value.
}


\references{Kratzer, G., Furrer, R., 2018. Information-Theoretic Scoring Rules to Learn Additive Bayesian Network Applied to Epidemiology. Preprint; Arxiv: stat.ML/1808.01126. 

  Further information about \bold{abn} can be found at:\cr
  \url{http://www.r-bayesian-networks.org}}

\author{Gilles Kratzer}

\examples{
## use built-in simulated data set

mydat <- ex0.dag.data[,c("b1","b2","b3","g1","b4","p2","p4")];## take a subset of cols

## setup distribution list for each node
mydists <- list(b1="binomial",
              b2="binomial",
              b3="binomial",
              g1="gaussian",
              b4="binomial",
              p2="poisson",
              p4="poisson"
             );
## null model - all independent variables
mydag.empty<-matrix(data=c(
                     0,0,0,0,0,0,0, # 
                     0,0,0,0,0,0,0, #
                     0,0,0,0,0,0,0, # 
                     0,0,0,0,0,0,0, # 
                     0,0,0,0,0,0,0, #
                     0,0,0,0,0,0,0, #
                     0,0,0,0,0,0,0  #
                     ), byrow=TRUE,ncol=7);
colnames(mydag.empty)<-rownames(mydag.empty)<-names(mydat);
## now repeat but include some dependencies
mydag<-mydag.empty;
mydag["b1","b2"]<-1; # b1<-b2 
mydag["b2","p4"]<-1; # b2<-p4
mydag["b2","g1"]<-1; # b2<-g1
mydag["g1","p2"]<-1; # g1<-p2
mydag["b3","g1"]<-1; # b3<-g1
mydag["b4","b1"]<-1; # b4<-b1
mydag["p4","g1"]<-1; # p4<-g1

## now fit the model to calculate its goodness of fit
myres.c <- fitabn.mle(dag.m=mydag,data.df=mydat,data.dists=mydists,centre=TRUE);

print(myres.c$mlik);

plotabn(dag.m = mydag,data.dists = mydists,fitted.values.abn.mle = myres.c$modes)
}

\keyword{models}
