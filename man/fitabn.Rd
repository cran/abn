% abn-internal.Rd --- 
% Author           : Fraser Lewis
% Last modified on : 06/08/2012 

\name{fitabn}
\alias{fitabn}
\encoding{latin1}
%- Also NEED an `\alias' for EACH other topic documented here.

\title{Fit an additive Bayesian network model}

\description{Fits an additive Bayesian network to observed data and is equivalent to Bayesian multi-dimensional regression modeling. Two numerical options are available, standard Laplace approximation or else an integrated nested Laplace approximation provided via a call to the R INLA library (see www.r-inla.org - this is not hosted on CRAN).}

\usage{

fitabn(dag.m=NULL, data.df=NULL, data.dists=NULL, group.var=NULL,cor.vars=NULL,
       create.graph=FALSE,compute.fixed=FALSE,mean=0, prec=0.001,loggam.shape=1,
       loggam.inv.scale=5e-05,verbose=FALSE, centre=TRUE,max.mode.error=10,
       max.iters=100,epsabs=1e-7,error.verbose=FALSE,epsabs.inner=1e-6,
       max.iters.inner=100,finite.step.size=1E-07,hessian.params=c(1E-04,1E-02),
       max.iters.hessian=10,max.hessian.error=1E-04,factor.brent=1E+02, 
       maxiters.hessian.brent=10,num.intervals.brent=100,min.pdf=1E-03,
       n.grid=100,std.area=TRUE, marginal.quantiles=c(0.025,0.25,0.5,0.75,0.975),
       max.grid.iter=1000,marginal.node=NULL, marginal.param=NULL,variate.vec=NULL)

}
  
%- maybe also `usage' for other objects documented here.
\arguments{ 
  \item{dag.m}{a matrix defining the network structure, a directed acyclic graph (DAG), see details for format. Note that colnames and rownames must be set}
  \item{data.df}{a data frame containing the data used for learning the network, binary variables must be declared as factors and no missing values all allowed in any variable}
  \item{data.dists}{a named list giving the distribution for each node in the network, see details}
  \item{group.var}{only applicable for mixed models and gives the column name in \code{data.df} of the grouping variable (which must be a factor denoting group membership)}
  \item{cor.vars}{a character vector giving the column names in data.df for which a mixed model should be used}
  \item{create.graph}{create an R graph class object which enables easy plotting of dag.m using \code{plot()}, requires Rgraphviz (hosted on bioconductor rather than CRAN).}
  \item{compute.fixed}{a logical flag, set to TRUE for computation of marginal posterior distributions, see details.}
  \item{mean}{the prior mean for all the Gaussian additive terms for each node.}
  \item{prec}{the prior precision for all the Gaussian additive terms for each node.}
  \item{loggam.shape}{the shape parameter in the Gamma distributed prior for the precision in any Gaussian nodes, also used for group level precision is applicable.}
  \item{loggam.inv.scale}{the inverse scale parameter in the Gamma distributed prior for the precision in any Gaussian nodes, also used for group level precision is applicable. }
  \item{verbose}{if TRUE then provides some additional output, in particular the code used to call INLA, if applicable.}
  \item{centre}{should the observations in each Gaussian node first be standarised to mean zero and standard deviation one}
  \item{max.mode.error}{if the estimated modes from INLA differ by a factor of max.mode.error or more from those computed internally then results from INLA are replaced by those computed internally. To force INLA to always be used then max.mode.error=100, to force INLA to never be used max.mod.error=0. See details.}  
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
\item{min.pdf}{the value of the posterior density function below which we stop the estimation, only used when computing marginals, see details.} 
\item{n.grid}{recompute density on an equally spaced grid with \code{n.grid} points.}
\item{std.area}{logical, should the area under the estimated posterior density be standardised to exactly one, useful for error checking.} 
\item{marginal.quantiles}{vector giving quantiles at which to compute the posterior marginal distribution at.}
\item{max.grid.iter}{gives number of grid points to estimate posterior density at when not explicitly specifying a grid, used to avoid excessively long computation.}
\item{marginal.node}{used in conjuction with \code{marginal.param} to allow bespoke estimate of a marginal density over a specific grid. value from 1 to number of nodes.} 
 \item{marginal.param}{used in conjunction with \code{marginal.node}. value of 1 is for intercept, see modes entry in results for appropriate number.}
\item{variate.vec}{a vector containing the places to evaluate the posterior marginal density, must be supplied if \code{marginal.node} is not null}
}

\details{
 The procedure \code{fitabn} fits an additive Bayesian network model to data where each node (variable - a column in data.df) can be either: presence/absence (Bernoulli); continuous (Gaussian); or an unbounded count (Poisson). The model comprises of a set of conditionally independent generalized linear regressions with or without random effects. Internal code is used by default for numerical estimation in nodes without random effects, and INLA is the default for nodes with random effects. This default behaviour can be over-ridden using max.mode.error. The default is \code{max.mode.error=10} which means that the modes estimated from INLA output must be within 10\% of those estimated using internal code, otherwise internal code is used rather than INLA. To force the use of INLA on all nodes use max.mode.error=100 which then ignores this check, to force the use of internal code then use max.mode.error=0. For the numerical reliability and perform of abn see \url{http://www.r-bayesian-networks.org/quality-assurance}. Generally speaking INLA can be extremely fast and accurate, but in a number of cases it can perform very poorly and so some care is required (which is why there is an internal check on the modes). Binary variables must be declared as factors with two levels, and the argument data.dists must be a list with named arguments, one for each of the variables in \code{data.df} (except a grouping variable - if present), where each entry is either "poisson","binomial", or "gaussian", see examples below. The "poisson" and "binomial" distributions use log and logit link functions respectively. Note that "binomial" here actually means only binary, one bernoulli trial per row in data.df. 

If the data are grouped into correlated blocks - where in a standard regression context a mixed model might be used - then a network comprising of one or more nodes where a generalized linear mixed model is used (but limited to only a single random effect). This is achieved by specifying parameters \code{group.var} and \code{cor.vars}, where the former defines the group membership variable which should be a factor indicating which observations belong to the same grouping. The parameter \code{cor.vars} is a character vector which contains the names of the nodes for which a mixed model should be used. For example, in some problems it may be appropriate for all variables (except group.var) in data.df to be parameterised as a mixed model while in others it may only be a single variable for which grouping adjustment is required (as the remainder of variables are covariates measured at group level).   

In the network structure definition, \code{dag.m}, each row represents a node in the network, and the columns in each row define the parents for that particular node, see the example below for the specific format.

If compute.fixed=TRUE then marginal posterior distributions for all parameters are computed. Note the current algorithm used to determine the evaluation grid is rather crude and may need to be manually refined using variate.vec (one parameter at a time) for publication quality density estimates. Note that a manual grid can only be used with internal code and not INLA (which uses its own grid). The end points are defined as where the value of the marginal density drops below a given threshold \code{pdf.min}. 

If create.graph=TRUE then the model definition matrix in \code{dag.m} is used to create an R graph object (of type graphAM-class). See ?"graph-class" for details and the Rgraphviz documentation (which is extensive). The main purpose of this is simply to allow easy visualisation of the DAG structure via the graphviz library. A graph plot can easily be created by calling the method plot() on this object (see example below). Note, however, that the fonts and choice of scaling used here may be far less visually optimal than using graphviz direct (e.g via \code{\link{tographviz}}) for publication quality output. Also, re-scaling the plotting window may not result in a callback to re-optimise the visual position of the nodes and edges, and so if the window is re-sized then re-run the plot command to re-draw to the new scale.   

When estimating the log marginal likelihood in models with random effects (using internal code rather than INLA) an attempt is made is minimize the error by comparing the estimates given between a 3pt and 5pt rule when estimating the Hessian in the Laplace approximation. The modes used in each case are identical. The first derivatives are computed using gsl's adaptive finite difference function and this is embedding inside the standard 3pt and 5pt rules for the second derivatives. In all cases a central difference approximation is tried first with a forward difference being a fall back (as the precision parameters are strictly positive). The error is minimised through choosing an optimal stepsize using gsl's Nelder-Mead optimization, and if this fails, (e.g. is greater than \code{max.hessian.error}) then the Brent-Dekker root bracketing method is used as a fall back. If the error cannot be reduced to below \code{max.hessian.error} then the stepsize which gave the lowest error during the searches (across potentially many different initial bracket choices) is used for the final Hessian evaluations in the Laplace approximation.

}

\value{A named list. One entry for each of the variables in \code{data.df} (excluding the grouping variable, if present) which contains an estimate of the log marginal likelihood for that individual node. An entry "mlik" which is the total log marginal likelihood for the full ABN model. A vector of error.codes - non-zero if a numerical error or warning occured, and a vector error.code.desc giving a text description of the error. A list "modes", which contains all the mode estimates for each parameter at each node. A vector called hessian accuracy which is the estimated local error in the log marginal likelihood for each node.  If compute.fixed=TRUE then a list entry called "marginals" which contains a named entry for every parameter in the ABN and each entry in this list is a two column matrix where the first column is the value of the marginal parameter, say x, and the second column is the respective density value, pdf(x). In addition a list called "marginal.quantiles" is produced giving the quantiles for each marginal posterior distribution. If create.graph=TRUE then an additional entry "graph" which is of type class graphAM-class is created. 

}


\references{Lewis FI, McCormick BJJ (2012). Revealing the complexity of health determinants in resource poor settings. \emph{American Journal Of Epidemiology}. DOI:10.1093/aje/KWS183).  

  Further information about \bold{abn} can be found at:\cr
  \url{http://www.r-bayesian-networks.org}}

\author{
  Fraser Lewis \email{fraseriain.lewis@uzh.ch}
}

%\seealso{\code{"\link[=fitabn.inla]{fitbn.inla}"}}

\examples{
\dontrun{
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

## now fit the model to calculate its goodness of fit
myres.c<-fitabn(dag.m=mydag.empty,data.df=mydat,data.dists=mydists,centre=TRUE,
                compute.fixed=FALSE);

print(myres.c$mlik);## log marginal likelihood goodness of fit for complete DAG

## now repeat but include some dependencies
mydag<-mydag.empty;
mydag["b1","b2"]<-1; # b1<-b2 
mydag["b2","p4"]<-1; # b2<-p4
mydag["b2","g1"]<-1; # b2<-g1
mydag["g1","p2"]<-1; # g1<-p2
mydag["b3","g1"]<-1; # b3<-g1
mydag["b4","b1"]<-1; # b4<-b1
mydag["p4","g1"]<-1; # p4<-g1

myres.c<-fitabn(dag.m=mydag,data.df=mydat,data.dists=mydists,centre=TRUE,
                compute.fixed=FALSE);

print(myres.c$mlik);## a much poorer fit than full independence DAG

## now also plot the graph via graphviz 
myres.c<-fitabn(dag.m=mydag,data.df=mydat,data.dists=mydists,centre=TRUE,
                create.graph=TRUE,compute.fixed=FALSE);

## ok for quick visualisation - but not publication quality 
plot(myres.c$graph);# see ?plot.graph for details.

## for publication quality may be better to use graphviz direct 
tographviz(dag.m=mydag,data.df=mydat,data.dists=mydists,outfile="graph.dot");
## and then process using graphviz tools e.g. on linux
system("dot -Tpdf -o graph.pdf graph.dot")
system("evince graph.pdf");

## a simple plot of some posterior densities
## the algorithm which chooses density points is very simple any may be rather sparse so also
## recompute the density over an equally spaced grid of 100 points between the two end points
## which had at f=min.pdf
myres.c<-fitabn(dag.m=mydag,data.df=mydat,data.dists=mydists,centre=TRUE,
                compute.fixed=TRUE,n.grid=100);

## repeat but use INLA for the numerics using max.mode.error=100
## as using internal code is the default here rather than INLA
myres.inla<-fitabn(dag.m=mydag,data.df=mydat,data.dists=mydists,centre=TRUE,
                compute.fixed=TRUE,max.mode.error=100);

print(names(myres.c$marginals));## gives all the different parameter names
## plot posterior densities
par(mfrow=c(2,2));
plot(myres.c$marginals$b1[["b1|(Intercept)"]],type="b",xlab="b1|(Intercept)",
main="Posterior Marginal Density",pch="+");
points(myres.inla$marginals$b1[["b1|(Intercept)"]],col="blue");
plot(myres.c$marginals$b2[["b2|p4"]],type="b",xlab="b2|p4",
main="Posterior Marginal Density",pch="+");
points(myres.inla$marginals$b2[["b2|p4"]],col="blue");
plot(myres.c$marginals$g1[["g1|precision"]],type="b",xlab="g1|precision",
main="Posterior Marginal Density",pch="+");
points(myres.inla$marginals$g1[["g1|precision"]],col="blue");
plot(myres.c$marginals$b4[["b4|b1"]],type="b",xlab="b4|b1",
main="Posterior Marginal Density",pch="+");
points(myres.inla$marginals$b4[["b4|b1"]],col="blue");

### A very simple mixed model example using built-in data

## specify dag - only two variables using subset of variables from ex3.dag.data
## both variables are assumed to need (separate) adjustment for the group variable
## i.e. a binomial glmm at each node

## model where b1<-b2
mydag<-matrix(data=c(
                     0,1, # b1
                     0,0  # b2
                     ), byrow=TRUE,ncol=2);

colnames(mydag)<-rownames(mydag)<-names(ex3.dag.data[,c(1,2)]);
## specific distributions
mydists<-list(b1="binomial",
              b2="binomial"
             );

## just compute marginal likelihood - use internal code via max.mode.error=0
## as using INLA is the default here.
myres.c<-fitabn(dag.m=mydag,data.df=ex3.dag.data[,c(1,2,14)],data.dists=mydists,
                group.var="group",cor.vars=c("b1","b2"),
                centre=TRUE,compute.fixed=FALSE,max.mode.error=0);
print(myres.c);## show all the output 

## compare modes for node b1 with glmer()
library(lme4);
m1<-glmer(b1~1+b2+(1|group),family="binomial",data=ex3.dag.data[,c(1,2,14)])
print(slot(m1,"fixef"));
print(1/unlist(VarCorr(m1)));
print(myres.c$modes$b1);## almost identical to lme4 n.b. glmer() gives variance
##                         fitabn gives precision=1/variance

## compare with INLA estimate
myres.inla<-fitabn(dag.m=mydag,data.df=ex3.dag.data[,c(1,2,14)],
                        data.dists=mydists,group.var="group",cor.vars=c("b1","b2"),
                        centre=TRUE,compute.fixed=FALSE,max.mode.error=100); 

## compare log marginal likelihoods for each node and total DAG - should be very similar
cbind(unlist(myres.inla[1:3]),unlist(myres.c[1:3]));

## now for marginals - INLA is strongly preferable for estimating marginals for nodes 
## with random effects as it is far faster, but may not be reliable
## see www.r-bayesian-networks.org/quality-assurance

# INLA's estimates of the marginals, using high n.grid=500
# as this makes the plots smoother - see below.
myres.inla<-fitabn(dag.m=mydag,data.df=ex3.dag.data[,c(1,2,14)],data.dists=mydists,
                        group.var="group",cor.vars=c("b1","b2"),
                        compute.fixed=TRUE,max.mode.error=100,n.grid=500);

## this is NOT recommended - marginal density estimation using fitabn in mixed models
## is really just for diagnostic purposes, better to use fitabn.inla() here
## but here goes...be patient                                                                                        
myres.c<-fitabn(dag.m=mydag,data.df=ex3.dag.data[,c(1,2,14)],data.dists=mydists,
                group.var="group",cor.vars=c("b1","b2"),compute.fixed=TRUE,max.mode.error=0);

## compare marginals between internal and INLA.   
par(mfrow=c(2,3))
## 5 parameters - two intercepts, one slope, two group level precisions
plot(myres.inla$marginals$b1[[1]],type="l",col="blue",lwd=1,pch="+");
points(myres.c$marginals$b1[[1]],col="brown",lwd=2);
plot(myres.inla$marginals$b1[[2]],type="l",col="blue",lwd=1,pch="+");
points(myres.c$marginals$b1[[2]],col="brown",lwd=2);
## the precision of group-level random effects
plot(myres.inla$marginals$b1[[3]],type="l",col="blue",xlim=c(0,2),lwd=1,pch="+");
points(myres.c$marginals$b1[[3]],col="brown",lwd=2);
plot(myres.inla$marginals$b2[[1]],type="l",col="blue",lwd=1,pch="+");
points(myres.c$marginals$b2[[1]],col="brown",lwd=2);
plot(myres.inla$marginals$b2[[1]],type="l",col="blue",lwd=1,pch="+");
points(myres.c$marginals$b2[[1]],col="brown",lwd=2);
## the precision of group-level random effects
plot(myres.inla$marginals$b2[[2]],type="l",col="blue",xlim=c(0,2),lwd=1,pch="+");
points(myres.c$marginals$b2[[2]],col="brown",lwd=2);

### these are very similar although not exactly identical

## use internal code but only to compute a single parameter over a specified grid
## this can be necessary if the simple auto grid finding functions does a poor job

myres.c<-fitabn(dag.m=mydag,data.df=ex3.dag.data[,c(1,2,14)],data.dists=mydists,
                group.var="group",
                cor.vars=c("b1","b2"),centre=FALSE,compute.fixed=TRUE,
                marginal.node=1,marginal.param=3,## precision term in node 1
                variate.vec=seq(0.05,1.5,len=25));

par(mfrow=c(1,2));
plot(myres.c$marginals[[1]],type="b",col="blue");## still fairly sparse
## an easy way is to use spline to fill in the density without recomputing other points
## provided the original grid is not too sparse.
plot(spline(myres.c$marginals[[1]],n=100),type="b",col="brown");
}

}

\keyword{models}
