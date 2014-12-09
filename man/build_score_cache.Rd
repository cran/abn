% buildscorecache.Rd --- 
% Author           : Fraser Lewis
% Last modified on : 03/08/2012 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
\name{buildscorecache}

\alias{buildscorecache}

\title{Build a cache of goodness of fit metrics for each node in a DAG, possibly subject to user defined restrictions}

\description{
  Iterates over all valid parent combinations - subject to ban, retain and max. parent limits - for each node, or a subset of nodes, and computes a cache of log marginal likelihoods. This cache can then be used in different DAG structural search algorithms. 
}

\usage{

buildscorecache(data.df=NULL, data.dists=NULL, group.var=NULL,cor.vars=NULL,
                dag.banned=NULL, dag.retained=NULL,max.parents=NULL,
                which.nodes=NULL,defn.res=NULL,dry.run=FALSE,
                max.mode.error=10,
                verbose=FALSE,centre=TRUE,mean=0, prec=0.001,loggam.shape=1,
                loggam.inv.scale=5e-05, max.iters=100,
                epsabs=1e-7,error.verbose=FALSE,
                epsabs.inner=1e-6,max.iters.inner=100,
                finite.step.size=1e-7,
                hessian.params=c(1E-04,1E-02),max.iters.hessian=10,
                max.hessian.error=5E-01,factor.brent=1E+02, 
                maxiters.hessian.brent=100,num.intervals.brent=100
)
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
 \item{max.mode.error}{if the estimated modes from INLA differ by a factor of max.mode.error or more from those computed internally then results from INLA are replaced by those computed internally. To force INLA to always be used then max.mode.error=100, to force INLA to never be used max.mod.error=0. See details.}  
  \item{verbose}{if TRUE then provides some additional output}
  \item{centre}{should the observations in each Gaussian node first be standarised to mean zero and standard deviation one, defaults to TRUE}
  \item{mean}{the prior mean for all the Gaussian additive terms for each node}
  \item{prec}{the prior precision for all the Gaussian additive term for each node}
  \item{loggam.shape}{the shape parameter in the Gamma distribution prior for the precision in a Gaussian node}
  \item{loggam.inv.scale}{the inverse scale parameter in the Gamma distribution prior for the precision in a Gaussian node}
   \item{max.iters}{total number of iterations allowed when estimating the modes in Laplace approximation}
  \item{epsabs}{absolute error when estimating the modes in Laplace approximation for models with no random effects.}
  \item{error.verbose}{logical,additional output in the case of errors occuring in the optimisation}
   \item{epsabs.inner}{absolute error in the maximisation step in the (nested) Laplace approximation for each random effect term}
   \item{max.iters.inner}{total number of iterations in the maximisation step in the nested Laplace approximation}
  \item{finite.step.size}{suggested step length used in finite difference estimation of the derivatives for the (outer) Laplace approximation when estimating modes}
  \item{hessian.params}{a numeric vector giving parameters for the adaptive algorithm which determines the optimal stepsize in the finite difference estimation of the hessian. First entry is the inital guess, second entry absolute error}
\item{max.iters.hessian}{integer, maximum number of iterations to use when determining an optimal finite difference approximation (nelder-mead)}
\item{max.hessian.error}{if the estimated log marginal likelihood when using an adaptive 5pt finite difference rule for the Hessian differs by more than max.hessian.error from when using an adaptive 3pt rule then continue to minimize the local error by switching to the Brent-Dekker root bracketing method, see details} 
\item{factor.brent}{if using Brent-Dekker root bracketing method then define the outer most interval end points as the best estimate of h (stepsize) from the nelder-mead as (h/factor.brent,h*factor.brent)}
\item{maxiters.hessian.brent}{maximum number of interations allowed in the Brent-Dekker method}
\item{num.intervals.brent}{the number of initial different bracket segments to try in the Brent-Dekker method}
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
This function is used to calculate all individual node scores (log marginal likelihoods). This cache can then be fed into a model search algorithm. This function is very similar to \code{\link{fitabn}} - see that help page for details of the type of models used and in particular data.dists specification - but rather than fit a single complete DAG \code{buildscorecache} iterates over all different parent combinations for each node.  There are three ways to customise the parent combinations through giving a matrix which contains arcs which are not allowed (banned), a matrix which contains arcs which must always be included (retained) and also a general complexity limit which restricts the maximum number of arcs allowed to terminate at a node (its number of parents), where this can differ from node to node. In these matrices, \code{dag.banned} and \code{dag.retained}, each row represents a node in the network, and the columns in each row define the parents for that particular node, see the example below for the specific format. If these are not supplied they are assumed to be empty matrices, i.e. no arcs banned or retained. Note that only rudimentary consistency checking is done here and some care should be taken not to provide conflicting restrictions in the ban and retain matrices.

The variable \code{which.nodes} is to allow the computation to be separated by node, for example over different cpus using say R CMD BATCH. This may useful and indeed likely essential with larger problems or those with random effects. Note that in this case the results must then be combined back into a list of identical format to that produced by an individual call to buildscorecache, comprising of all nodes (in same order as the columns in data.df) before sending to any search routines. Using \code{dry.run} can be useful here.  

The numerical routines used here are identical to those in \code{\link{fitabn}} and see that help page for further details and also the quality assurance section on the \href{www.r-bayesian-networks.org}{abn website} for more details.
 
}

\references{Lewis FI, McCormick BJJ (2012). Revealing the complexity of health determinants in resource poor settings. \emph{American Journal Of Epidemiology}. DOI:10.1093/aje/KWS183).  

  Further information about \bold{abn} can be found at:\cr
  \url{http://www.r-bayesian-networks.org}}

\author{Fraser Ian Lewis}

\examples{
\dontrun{
#################################################################
## Example 1
#################################################################
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
                     dag.banned=ban, dag.retained=retain,max.parents=max.par
                     );

## repeat but using R-INLA. The mlik's should be virtually identical.
## now build cache
res.inla<-buildscorecache(data.df=mydat,data.dists=mydists,
                     dag.banned=ban, dag.retained=retain,max.parents=max.par,
                     max.mode.error=100);

## plot comparison - very similar
plot(res.c$mlik,res.inla$mlik,pch="+");abline(0,1);


#################################################################
## Example 2 - much bigger problem using glm - make take a while
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
max.par<-list("b1"=4,"g1"=4,"p1"=4,"b2"=4,"g2"=4,"p2"=4,"b3"=4,
              "g3"=4,"p3"=4,"b4"=4,"g4"=4,
               "p4"=4,"b5"=4,"g5"=4,"p5"=4,"b6"=4,"g6"=4,"p6"=4);

## no explicit ban or retain restrictions set so dont need to supply ban
##  or retain matrices

## now build cache using internal code just for nodes 1,2 and 3
## e.g. "b1", "p1" and "g1" 
mycache.c<-buildscorecache(data.df=mydat,data.dists=mydists,
                         max.parents=max.par, which.nodes=c(1:3));

###################################################################
## Example 3 - grouped data - random effects example e.g. glmm
###################################################################

mydat<-ex3.dag.data;## this data comes with abn see ?ex3.dag.data

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

## in this example INLA is used as default since these are glmm nodes
## when running this at node-parent combination 71 the default accuracy check on the 
## INLA modes is exceeded (default is a max. of 10 percent difference from
## modes estimated using internal code) and a message is given that internal code
## will be used in place of INLA's results.

mycache<-buildscorecache(data.df=mydat,data.dists=mydists,group.var="group",
                         cor.vars=c("b1","b2","b3","b4","b5","b6","b7",
                                    "b8","b9","b10",
                                    "b11","b12","b13"),
                         max.parents=max.par, which.nodes=c(1));


}

}

\keyword{laplace}
