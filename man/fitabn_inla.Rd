% abn-internal.Rd --- 
% Author           : Fraser Lewis
% Last modified on : 06/08/2012 

\name{fitabn.inla}
\alias{"fitabn.inla"}
\encoding{latin1}
%- Also NEED an `\alias' for EACH other topic documented here.
\alias{fitabn.inla}
\title{Fit an additive Bayesian network model using INLA}

\description{Fits an additive Bayesian network to observed data and is equivalent to Bayesian multi-dimensional regression modeling. R-INLA is used for the model fitting and must be available.}

\usage{

fitabn.inla(dag.m=NULL, data.df=NULL, data.dists=NULL, group.var=NULL,cor.vars=NULL,
            create.graph=FALSE,compute.fixed=FALSE,mean=0, prec=0.001,
            loggam.shape=1,loggam.inv.scale=5e-05, verbose=FALSE, centre=TRUE)

}

  
%- maybe also `usage' for other objects documented here.
\arguments{ 
  \item{dag.m}{a matrix defining the network structure, a directed acyclic graph (DAG), see details for format. Note that colnames and rownames must be set}
  \item{data.df}{a data frame containing the data used for learning the network, binary variables must be declared as factors and no missing values}
  \item{data.dists}{a named list giving the distribution for each node in the network, see details}
  \item{group.var}{only applicable for mixed models and gives the column name in \code{data.df} of the grouping variable (which must be a factor denoting group membership)}
  \item{cor.vars}{a character vector giving the column names in data.df for which a mixed model should be used to adjust for within group correlation}
  \item{create.graph}{create an R graph class object which enables easy plotting of dag.m using \code{plot()}, requires Rgraphviz.}
  \item{compute.fixed}{a logical flag, set to TRUE for computation of marginal posterior distributions as well as the (log) marginal likelihood.}
  \item{mean}{the prior mean of the Gaussian additive terms for each node}
  \item{prec}{the prior precision of the Gaussian additive term for each node}
  \item{loggam.shape}{the shape parameter in the Gamma distributed prior for the precision in any Gaussian nodes}
  \item{loggam.inv.scale}{the inverse scale parameter in the Gamma distributed prior for the precision in any Gaussian nodes}
  \item{verbose}{if true then provides some additional output}
  \item{centre}{should the observations in each Gaussian node first be standarised to mean zero and standard deviation one}
  }

\details{
 The procedure \code{fitabn.inla} fits an additive Bayesian network model to data where each node (variable - a column in data.df) can be either: presence/absence (Bernoulli); continuous (Gaussian); or an unbounded count (Poisson). The model comprises of a set of conditionally independent generalized linear regressions. R-INLA is used for the model fitting - this function is effectively just a wrapper to \code{inla()}. If the INLA library is not available this function with give an error. R-INLA can be downloaded from \url{http://www.r-inla.org/downloadavailable}. Currently implemented are binary (binomial data must be presented as individual bernoulli trials in each row of data.df), Gaussian or Poisson likelihoods. Binary variables must be declared as factors with two levels, and the argument data.dists must be a list with named arguments, one for each of the variables in \code{data.df} (except a grouping variable - if present), where each entry is either "poisson","binomial", or "gaussian", see examples below. The "poisson" and "binomial" distributions use log and logit link functions respectively. 

If the data is grouped into correlated blocks - where in a standard regression context a mixed model might be used - then a network comprising of one or more nodes where a generalized linear mixed model is used. This is achieved by specifying parameters \code{group.var} and \code{cor.vars}, where the former defines the group membership variable which should be a factor indicating which observations below to the same grouping. The parameter \code{cor.vars} is a character vector which contains the names of the nodes for which a mixed model should be used, for example in some problems is may be appropriate for all variables (except group.var) in data.df to be parameterised as a mixed model which in others it may only be a single variable for which grouping adjustment is required.   

In the network structure definition, \code{dag.m}, each row represents a node in the network, and the columns in each row define the parents for that particular node, see the example below for the specific format.

If compute.fixed=TRUE then marginal posterior distributions for all parameters are computed. Note that for models without random effects R-INLA can be very slow in such cases and may be less reliable than using \code{\link{fitabn}}.   

If create.graph=TRUE then the model definition matrix in \code{dag.m} is used to create an R graph object (of type graphAM-class). See ?"graph-class" for details and the Rgraphviz documentation (which is extensive). The main purpose of this is simply to allow easy visualisation of the DAG structure via the graphviz library. A graph plot can easily be created by calling the method plot() on this object (see example below). Note, however, that the fonts and choice of scaling used here may be far less visually optimal than using graphviz direct (e.g via \code{\link{tographviz}}) for publication quality output. Also, re-scaling the plotting window does not result in a callback to re-optimise the visual position of the nodes and edges, and so if the window is re-sized then re-run the plot command to re-draw to the new scale.   

}

\value{A named list. One entry for each of the variables in \code{data.df} (excluding the grouping variable, if present) which contains an estimate of the log marginal likelihood for that individual node. An entry "mlik" which is the total log marginal likelihood for the full DAG model.  If compute.fixed=TRUE then a list entry called "marginals" which contains a named entry for every parameter in the DAG and each entry in this list is a two column matrix where the first column is the value of the marginal parameter, say x, and the second column is the respective density value, pdf(x). If create.graph=TRUE then an additional entry "graph" which is of type class graphAM-class is created.  

}


\references{ Rue H., Martino S. and Chopin N. (2009): Approximate Bayesian Inference for Latent Gaussian Models Using Integrated Nested Laplace Approximations (with discussion). Journal of the Royal Statistical Society, Series B, 71, 319-392

  Further information about \bold{abn} can be found at:\cr
  \url{http://www.r-bayesian-networks.org}}

\author{
  Fraser Lewis \email{fraseriain.lewis@uzh.ch}
}

\seealso{\code{"\link[=fitabn]{fitbn}"}}

\examples{

\dontrun{## 
## use built-in simulated data set

mydat<-ex0.dag.data[,c("b1","b2","b3","g1","b4","p2","p4")];## take a subset of cols

## setup distribution list for each node
mydists<-list(b1="binomial",
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

## fit model to get goodness of fit
myres.inla<-fitabn.inla(dag.m=mydag.empty,data.df=mydat,data.dists=mydists,centre=TRUE,
                compute.fixed=FALSE);

print(myres.inla$mlik);## log marginal likelihood goodness of fit for complete DAG

## now repeat but include some dependencies
mydag<-mydag.empty;
mydag["b1","b2"]<-1; # b1<-b2 
mydag["b2","p4"]<-1; # b2<-p4
mydag["b2","g1"]<-1; # b2<-g1
mydag["g1","p2"]<-1; # g1<-p2
mydag["b3","g1"]<-1; # b3<-g1
mydag["b4","b1"]<-1; # b4<-b1
mydag["p4","g1"]<-1; # p4<-g1

myres.inla<-fitabn.inla(dag.m=mydag,data.df=mydat,data.dists=mydists,centre=TRUE,
                compute.fixed=FALSE);

print(myres.inla$mlik);## a much poorer fit than full independence DAG

## now also plot the graph via graphviz 
myres.inla<-fitabn.inla(dag.m=mydag,data.df=mydat,data.dists=mydists,centre=TRUE,
                create.graph=TRUE,compute.fixed=FALSE);

## ok for quick visualisation - but not publication quality 
plot(myres.inla$graph);# see ?plot.graph for details.

## for publication quality may be better to use graphviz direct 
tographviz(dag.m=mydag,data.df=mydat,data.dists=mydists,outfile="graph.dot");
## and then process using graphviz tools e.g. on linux
system("dot -Tpdf -o graph.pdf graph.dot")
system("evince graph.pdf");

## a simple plot of some posterior densities
myres.inla<-fitabn.inla(dag.m=mydag,data.df=mydat,data.dists=mydists,centre=TRUE,
                compute.fixed=TRUE,create.graph=TRUE);

## compare with non-INLA version
myres.c<-fitabn(dag.m=mydag,data.df=mydat,data.dists=mydists,centre=TRUE,
                compute.fixed=TRUE,n.grid=100);

print(names(myres.c$marginals));## gives all the different parameter names
## the plot below shows only red as the lines are virtually on top of each other
par(mfrow=c(2,2));
plot(myres.c$marginals[["b1|(Intercept)"]],type="b",xlab="b1|(Intercept)", 
main="Posterior Marginal Density",pch="+");
points(myres.inla$marginals[["b1|(Intercept)"]],col="red");
plot(myres.c$marginals[["b2|p4"]],type="b",xlab="b2|p4",
main="Posterior Marginal Density",pch="+");
points(myres.inla$marginals[["b2|p4"]],col="red");
plot(myres.c$marginals[["g1|precision"]],type="b",xlab="g1|precision",
main="Posterior Marginal Density",pch="+");
points(myres.inla$marginals[["g1|precision"]],col="red");
plot(myres.c$marginals[["b4|b1"]],type="b",xlab="b4|b1",
main="Posterior Marginal Density",pch="+");
points(myres.inla$marginals[["b4|b1"]],col="red");

### A very simple mixed model example using built-in data

## specify dag - only two variables using subset of variables from ex3.dag.data

mydag<-matrix(data=c(
                     0,1, # b1
                     0,0  # b2
                     ), byrow=TRUE,ncol=2);

colnames(mydag)<-rownames(mydag)<-names(ex3.dag.data[,c(1,2)]);
## specific distributions
mydists<-list(b1="binomial",
              b2="binomial"
             );

## the variable "group" in grouped.dag.data contains the group membership information
## (indices 1 through 50
## fit model to node 1 e.g. b1 as a binary mixed model with single covariate (e.g. glmm)
## for node 2 e.g. g1 just fit standard gaussian glm 
## just compute marginal likelihood
myres.inla<-fitabn.inla(dag.m=mydag,data.df=ex3.dag.data[,c(1,2,14)],data.dists=mydists,
                        group.var="group",cor.vars=c("b1","b2"),
                        centre=TRUE,compute.fixed=TRUE);

## what marginals do we have?
print(names(myres.inla$marginals));
plot(myres.inla$marginals[["b1|Groups-precision"]]);# marg. density for precision of random effects
## the x range here is too wide here but this determined by INLA and not the abn library

}}

\keyword{models}
