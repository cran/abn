% abn-internal.Rd --- 
% Author           : Fraser Lewis
% Last modified on : 03/08/2012 

\name{buildscorecache.inla}

\alias{buildscorecache.inla}

\title{Build a cache of goodness of fit metrics for each node in a DAG using R-INLA, possibly subject to user defined restrictions}

\description{
  Iterates over all valid parent combinations - subject to ban, retain and max. parent limits - for each node, or a subset of nodes, and computes a cache of log marginal likelihoods. This cache is then used in different DAG structural search algorithms. The R-INLA library is used for the numerics.
}

\usage{

buildscorecache.inla(data.df=NULL, data.dists=NULL, ntrials=NULL, exposure=NULL, 
                     group.var=NULL,cor.vars=NULL,dag.banned=NULL, 
                     dag.retained=NULL,max.parents=NULL,which.nodes=NULL,
                     defn.res=NULL,dry.run=FALSE,verbose=FALSE,centre=TRUE,
                     mean=0, prec=0.001,loggam.shape=1,loggam.inv.scale=5e-05)
}
\arguments{
  \item{data.df}{a data frame containing the data used for learning each node, binary variables must be declared as factors}
  \item{data.dists}{a named list giving the distribution for each node in the network, see details} 
  \item{ntrials}{a numeric vector giving the total number of trials, only applicable if the data comprise one or more binary variables. see details}
  \item{exposure}{a numeric vector, giving the unit of exposure, only applicable if the data comprise one or more poisson variables, see details}
  \item{group.var}{only applicable for nodes to be fitted as a mixed model and gives the column name in \code{data.df} of the grouping variable which must be a factor denoting group membership}
  \item{cor.vars}{a character vector giving the column names in data.df for which a mixed model should be used to adjust for within group correlation}
  \item{dag.banned}{a matrix defining which arcs are not permitted - banned - see details for format. Note that colnames and rownames must be set. }
  \item{dag.retained}{a matrix defining which arcs are must be retained in any model search, see details for format. Note that colnames and rownames must be set}
  \item{max.parents}{a constant or named list giving the maximum number of parents allowed, the list version allows this to vary per node.}
  \item{which.nodes}{a vector giving the column indices of the variables to be included, if ignored all variables are included}
  \item{defn.res}{an optional user-supplied list of child and parent combinations, see details}
  \item{dry.run}{if TRUE then a list of the child nodes and parent combinations are returned but without estimation of node scores (log marginal likelihoods)}
  \item{verbose}{if true then provides some additional output, in particular the call used to inla()}
  \item{centre}{should the observations in each Gaussian node first be standarised to mean zero and standard deviation one, defaults to TRUE}
  \item{mean}{the prior mean of the Gaussian additive terms for each node}
  \item{prec}{the prior precision of the Gaussian additive term for each node}
  \item{loggam.shape}{the shape parameter in the Gamma distribution prior for the precision in a Gaussian node}
  \item{loggam.inv.scale}{the inverse scale parameter in the Gamma distribution prior for the precision in a Gaussian node}
  }


\value{A named list. In addition to those members below this list may also contain a vector error.indexes if any of the node combinations could not be reliably estimated, see details.
\item{children}{a vector of the child node indexes (from 1) corresponding to the columns in data.df}
\item{node.defn}{a matrix giving the parent combination}
\item{mlik}{log marginal likelihood value for each node combination.}
\item{data.df}{a version of the original data (for internal use only in other functions such as mostprobable()).}
}

\details{
This function is used to calculate all individual node scores (log marginal likelihoods) using calls to R-INLA (this R library must be available - R-INLA can be downloaded from \url{http://www.r-inla.org/downloadavailable}). This cache can then be fed into a model search algorithm. This function is very similar to \code{\link{fitabn.inla}} - see that help page for details of the type of models used and in particular data.dists specification - but rather than fit a single complete DAG it iterates over all different parent combinations for each node.  There are three ways to customise the parent combinations through giving a matrix which contains arcs which are not allowed (banned), a matrix which contains arcs which must always be included (retained) and also a general complexity limit which restricts the maximum number of arcs allowed to terminate at a node (its number of parents), where this can differ from node to node. In these matrices, \code{dag.banned} and \code{dag.retained}, each row represents a node in the network, and the columns in each row define the parents for that particular node, see the example below for the specific format. If these are not supplied they are assumed to be empty matrices, i.e. no arcs banned or retained.

The variable \code{which.nodes} is to allow the computation to be separated by node, for example over different cpus using say R CMD BATCH. This may useful and indeed likely essential with larger problems. Note that the results must then be combined back into a list of identical format to that produced by an individual call to buildscorecache, comprising of all nodes (in same order as the columns in data.df) before sending to any search routines. The computation of the node cache is only a small faction of the overal computation required to identify high scoring full DAG models.

This function can provide a useful comparator to \code{\link{buildscorecache}} but may be less suitable for models without random effects at each node. It is certainly faster for models with random effects, possibly much faster, although does not provide any estimates of accuracy. Also note that in the (unlikely) event that INLA crashes due to a numerical error/difficulty then the whole run will terminate (which should not happen with \code{\link{buildscorecache}}). See the quality assurance section on the \href{www.r-bayesian-networks.org}{abn website} for more details of numerical comparisons.  
 
}

 \references{ Further information about \bold{abn} can be found at:\cr
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

res.c<-buildscorecache.inla(data.df=mydat,data.dists=mydists,
                     dag.banned=ban, dag.retained=retain,max.parents=max.par,
                     verbose=FALSE,centre=TRUE);

################################################################################################
## Example 2 glmm
################################################################################################

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


mycache.inla<-buildscorecache.inla(data.df=mydat,data.dists=mydists,group.var="group",
                         cor.vars=c("b1","b2","b3","b4","b5","b6","b7","b8","b9","b10","b11","b12","b13"),
                         max.parents=max.par, which.nodes=c(1),
                         verbose=FALSE,centre=TRUE);


}

}

\keyword{laplace}
