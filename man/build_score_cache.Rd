% abn-internal.Rd --- 
% Author           : Fraser Lewis
% Last modified on : 03/08/2012 

\name{buildscorecache}

\alias{buildscorecache}

\title{Build a cache of goodness of fit metrics for each node in a DAG, possibly subject to user defined restrictions}

\description{
  Iterates over all valid parent combinations - subject to ban, retain and max. parent limits - for each node, or a subset of nodes, and computes a cache of log marginal likelihoods. This cache is then used in different DAG structural search algorithms. 
}

\usage{

buildscorecache(data.df=NULL, data.dists=NULL, group.var=NULL,cor.vars=NULL,
                            dag.banned=NULL, dag.retained=NULL,max.parents=NULL,
                            which.nodes=NULL,defn.res=NULL,dry.run=FALSE,
                            verbose=FALSE,centre=TRUE,mean=0, prec=0.001,loggam.shape=1,
                            loggam.inv.scale=5e-05, max.iters=100,
                            epsabs=1e-7,error.verbose=FALSE,term.output.freq=100,
                            epsabs.inner=1e-6,max.iters.inner=100,
                            finite.step.size=1e-7,
                            hessian.params=c(1E-04,1E-02),max.iters.hessian=10)

}
\arguments{
  \item{data.df}{a data frame containing the data used for learning each node, binary variables must be declared as factors}
  \item{data.dists}{a named list giving the distribution for each node in the network, see details} 
  \item{group.var}{only applicable for nodes to be fitted as a mixed model and gives the column name in \code{data.df} of the grouping variable which must be a factor denoting group membership}
  \item{cor.vars}{a character vector giving the column names in \code{data.df} for which a mixed model should be used to adjust for within group correlation}
  \item{dag.banned}{a matrix defining which arcs are not permitted - banned - see details for format. Note that colnames and rownames must be set. }
  \item{dag.retained}{a matrix defining which arcs are must be retained in any model search, see details for format. Note that colnames and rownames must be set}
  \item{max.parents}{a constant or named list giving the maximum number of parents allowed, the list version allows this to vary per node.}
  \item{which.nodes}{a vector giving the column indices of the variables to be included, if ignored all variables are included}
  \item{defn.res}{an optional user-supplied list of child and parent combinations, see details}
  \item{dry.run}{if TRUE then a list of the child nodes and parent combinations are returned but without estimation of node scores (log marginal likelihoods)}
  \item{verbose}{if TRUE then provides some additional output}
  \item{centre}{should the observations in each Gaussian node first be standarised to mean zero and standard deviation one, defaults to TRUE}
  \item{mean}{the prior mean of the Gaussian additive terms for each node}
  \item{prec}{the prior precision of the Gaussian additive term for each node}
  \item{loggam.shape}{the shape parameter in the Gamma distribution prior for the precision in a Gaussian node}
  \item{loggam.inv.scale}{the inverse scale parameter in the Gamma distribution prior for the precision in a Gaussian node}
  \item{max.iters}{total number of iterations in the maximsation step in the Laplace approximation}
  \item{epsabs}{absolute error in the maximsation step in the Laplace approximation}
  \item{error.verbose}{additional output in the case of errors occuring in the optimisation}
  \item{term.output.freq}{a progress counter of how many parent combinations have so far been computed}
  \item{epsabs.inner}{absolute error in the maximisation step in the nested Laplace approximation for each random effect term}
  \item{max.iters.inner}{total number of iterations in the maximisation step in the nested Laplace approximation}
  \item{finite.step.size}{suggested step length used in finite difference estimation of the derivatives and hessian for the (outer) Laplace approximation - only applicable for nodes to be fitted as mixed models}
  \item{hessian.params}{a numeric vector giving parameters for the adaptive algorithm which determines the optimal stepsize in the finite difference estimation of the hessian. First entry is the inital guess, second entry absolute error, see details.} 
\item{max.iters.hessian}{integer, maximum number of iterations to use when determining an optimal finite difference approximation}
}


\value{A named list. 
\item{children}{a vector of the child node indexes (from 1) corresponding to the columns in data.df (ignoring any grouping variable)}
\item{node.defn}{a matrix giving the parent combination}
\item{mlik}{log marginal likelihood value for each node combination. If the model cannot be fitted then NA is returned. }
\item{error.code}{if non-zero then either the root finding algorithm (glm nodes) or the maximisation algorithm (glmm nodes) terminated in an unusual way suggesting a possible unreliable result, or else the finite difference hessian estimation produced and error or warning (glmm nodes)}
\item{error.code.desc}{a textual description of the error.code}
\item{hessian.accuracy}{An estimate of the error in the final mlik value for each parent combination - this is the absolute difference between two different adaptive finite difference rules where each computes the mlik value.}
\item{data.df}{a version of the original data (for internal use only in other functions such as mostprobable()).}
}

\details{
This function is used to calculate all individual node scores (log marginal likelihoods). This cache can then be fed into a model search algorithm. This function is very similar to \code{\link{fitabn}} - see that help page for details of the type of models used and in particular data.dists specification - but rather than fit a single complete DAG it iterates over all different parent combinations for each node.  There are three ways to customise the parent combinations through giving a matrix which contains arcs which are not allowed (banned), a matrix which contains arcs which must always be included (retained) and also a general complexity limit which restricts the maximum number of arcs allowed to terminate at a node (its number of parents), where this can differ from node to node. In these matrices, \code{dag.banned} and \code{dag.retained}, each row represents a node in the network, and the columns in each row define the parents for that particular node, see the example below for the specific format. If these are not supplied they are assumed to be empty matrices, i.e. no arcs banned or retained. Note that only rudimentary consistency checking is done here and some care should be taken not to provide conflicting restrictions in the ban and retain matrices.

The variable \code{which.nodes} is to allow the computation to be separated by node, for example over different cpus using say R CMD BATCH. This may useful and indeed likely essential with larger problems or those with random effects. Note that the results must then be combined back into a list of identical format to that produced by an individual call to buildscorecache, comprising of all nodes (in same order as the columns in data.df) before sending to any search routines. Using \code{dry.run} can be useful here.

Using this function as opposed to the similar \code{\link{buildscorecache.inla}} is recommended as this performs a number of validation checks on the marginal likelihood estimates, and should gracefully cope with models which cannot be fitted reliably without terminating part way through the cache computation. Generally speaking, if the node combination can be fitted using one of R's maximum likelihood routines (e.g. glm() or glmer()) then we would expect the numerical scheme used to implement to Laplace approximation to also work sensibly. If the internal Laplace approx fails then \code{\link{buildscorecache.inla}} may help in some cases as this uses a different numerical implementation. Empirical tests suggest that differences between the internal estimates and INLA estimates appear generally very small, although it is not clear when there are larger differences which is the more reliable. The inclusion of random effects makes the numerical problem vastly harder than when only fixed effects are present. The internal code provide some diagnostic information on the reliability of each estimate. In particular, as the hessian may be sensitive to the choice of finite difference step size two, different adaptive methods are used (different combinations of three and five point rules) and the absolute error in the Laplace approximation between these different approaches gives an local error estimate (hessian.accuracy) in the output. See the quality assurance section on the \href{www.r-bayesian-networks.org}{abn website} for more details. Note that the parameter finite.step.size is not the step size actually used in the finite difference computation, but rather only an initial suggestion provided to the GSL numerical differentiation functions which then determine an optimal step size. Secondly, the absolute error given in the (second) hessian.params argument is not the required minimum error in the mlik value, but rather the error in the optimal step size to be used in the Laplace approximation. Generally, the smaller this parameter then the smaller the Laplace error estimate but this not the criteria used in the minimisation and it may be impossible to reduce the mlik error this beyond a fixed minimum bound irrespective of stepsize. This can usually be seen from the results (in hessian.accuracy) as most nodes will likely have very lower mlik error but some may always have a larger error estimate irrespective of the numerical parameters used.    
 
}

\references{Lewis FI, McCormick BJJ (2012). Revealing the complexity of health determinants in resource poor settings. \emph{American Journal Of Epidemiology}. DOI:10.1093/aje/KWS183).  

  Further information about \bold{abn} can be found at:\cr
  \url{http://www.r-bayesian-networks.org}}
\author{
  Fraser Lewis \email{fraseriain.lewis@uzh.ch}.
}

\examples{
\dontrun{
## example 1
mydat<-ex0.dag.data[,c("b1","b2","g1","g2","b3","g3")];## take a subset of cols

## setup distribution list for each node
mydists<-list(b1="binomial",
              b2="binomial",
              g1="gaussian",
              g2="gaussian",
              b3="binomial",
              g3="gaussian"
             );

ban<-matrix(rep(0,dim(mydat)[2]^2),ncol=dim(mydat)[2]);# ban nothing
colnames(ban)<-rownames(ban)<-names(mydat); #names must be set
ban["b1","b2"]<-1; # now ban arc from b2 to b1 
retain<-matrix(rep(0,dim(mydat)[2]^2),ncol=dim(mydat)[2]);# retain nothing
colnames(retain)<-rownames(retain)<-names(mydat); #names must be set
retain["g1","g3"]<-1; # always retain arc from g3 to g1
# parent limits
max.par<-list("b1"=4,"b2"=4,"g1"=4,"g2"=0,"b3"=4,"g3"=4);

## now build cache of scores (goodness of fits for each node)

res.c<-buildscorecache(data.df=mydat,data.dists=mydists,
                     dag.banned=ban, dag.retained=retain,max.parents=max.par,
                     verbose=FALSE,centre=TRUE);

## repeat but using R-INLA. The mlik's should be virtually identical.
## now build cache
res.inla<-buildscorecache.inla(data.df=mydat,data.dists=mydists,
                     dag.banned=ban, dag.retained=retain,max.parents=max.par,
                     verbose=FALSE,centre=TRUE);

## plot comparison - very similar
plot(res.c$mlik,res.inla$mlik,pch="+");abline(0,1);


#################################################################
## example 2 - much bigger problem using glm - make take a while
#################################################################

mydat<-ex2.dag.data;## this data comes with abn see ?ex2.dag.data

## setup distribution list for each node
mydists<-list(b1="binomial",
              g1="gaussian",
              p1="poisson",
              b2="binomial",
              g2="gaussian",
              p2="poisson",
              b3="binomial",
              g3="gaussian",
              p3="poisson",
              b4="binomial",
              g4="gaussian",
              p4="poisson",
              b5="binomial",
              g5="gaussian",
              p5="poisson",
              b6="binomial",
              g6="gaussian",
              p6="poisson"
             );

## parent limits
max.par<-list("b1"=4,"g1"=4,"p1"=4,"b2"=4,"g2"=4,"p2"=4,"b3"=4,"g3"=4,"p3"=4,"b4"=4,"g4"=4,
"p4"=4,"b5"=4,"g5"=4,"p5"=4,"b6"=4,"g6"=4,"p6"=4);

## no explicit ban or retain restrictions set so dont need to supply ban or retain matrices

## now build cache using internal code just for nodes 1,2 and 3 e.g. "b1", "p1" and "g1" 
mycache.c<-buildscorecache(data.df=mydat,data.dists=mydists,
                         max.parents=max.par, which.nodes=c(1:3),
                         verbose=FALSE,centre=TRUE);

## smaller example to allow comparison with INLA 
mydat<-ex2.dag.data;## this data comes with abn see ?ex2.dag.data

## setup distribution list for each node
mydists<-list(b1="binomial",
              g1="gaussian",
              p1="poisson",
              b2="binomial",
              g2="gaussian",
              p2="poisson",
              b3="binomial",
              g3="gaussian",
              p3="poisson",
              b4="binomial",
              g4="gaussian",
              p4="poisson",
              b5="binomial",
              g5="gaussian",
              p5="poisson",
              b6="binomial",
              g6="gaussian",
              p6="poisson"
             );

max.par<-2;## dont use list use 2 parent max for all nodes

## no explicit ban or retain restrictions set so dont need to supply ban or retain matrices

## now build cache using internal code just for nodes 1,2 and 3 e.g. "b1", "p1" and "g1"
mycache.c<-buildscorecache(data.df=mydat,data.dists=mydists,
                         max.parents=max.par, which.nodes=c(1:3),
                         verbose=FALSE,centre=TRUE);
### This can take a LONG time
mycache.inla<-buildscorecache.inla(data.df=mydat,data.dists=mydists,
                         max.parents=max.par, which.nodes=c(1:3),
                         verbose=FALSE,centre=TRUE);
plot(mycache.inla$mlik,mycache.c$mlik,pch="+"); abline(0,1);

## all except 10 nodes-parent combinations are almost identical in agreement. A validation 
## of these can be found in the QA section on www.r-bayesian-networks.org. For this
## particular data in each case the value produced by the internal abn code is correct
## where as that from INLA is not. This particular data is ackward for INLA to deal with 
## due to a lack of latent (random effect). See exceptions on INLA website.
################################################################################################
## Example 3 - grouped data - random effects example e.g. glmm
################################################################################################

mydat<-ex3.dag.data;## this data comes with abn see ?ex1.dag.data

mydists<-list(b1="binomial",
              b2="binomial",
              b3="binomial",
              b4="binomial",
              b5="binomial",
              b6="binomial",
              b7="binomial",
              b8="binomial",
              b9="binomial",
              b10="binomial",
              b11="binomial",
              b12="binomial",
              b13="binomial"
             );
max.par<-2;

mycache.c<-buildscorecache(data.df=mydat,data.dists=mydists,group.var="group",
                         cor.vars=c("b1","b2","b3","b4","b5","b6","b7","b8","b9","b10",
                                    "b11","b12","b13"),
                         max.parents=max.par, which.nodes=c(1),finite.step.size=1e-07,
                         hessian.params=c(1E-04,1E-02),
                         term.output.freq=1,max.iters.hessian=10,
                         verbose=FALSE,centre=TRUE);

## compare with INLA. INLA is faster here
mycache.inla<-buildscorecache.inla(data.df=mydat,data.dists=mydists,group.var="group",
                         cor.vars=c("b1","b2","b3","b4","b5","b6","b7","b8","b9","b10",
                                    "b11","b12","b13"),
                         max.parents=max.par, which.nodes=c(1),
                         verbose=FALSE,centre=TRUE);

## a very close match - open question as to which is more accurate. Some small differences.
cbind(mycache.c$mlik,mycache.inla$mlik);

## The internal code version provides an estimate of the error in mlik values e.g. an idea
## of its robustness. These should be considered relative to the other node-parent 
## combinations. here these are mostly very small
plot(mycache.c$hessian.accuracy);## this is the absolute (estimated) error in the mlik value
## but as seems common some combinations have a rather larger error than the rest
plot(log10(mycache.c$hessian.accuracy));


}

}

\keyword{laplace}
