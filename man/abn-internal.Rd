% abn-internal.Rd --- 
% Author           : Fraser Lewis
% Last modified on : 26/09/2014, 30/05/2019 , 11/11/2019
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
\name{abn-internal}

\alias{calc.node.inla.glm}
\alias{calc.node.inla.glmm}
\alias{buildScoreCache.mle}
\alias{buildScoreCache.bayes}
\alias{fitAbn.bayes}
\alias{fitAbn.mle}

\title{abn internal functions}

\description{
These are functions for either internal use or error checking and are
not meant to be called directly by the user.
}

\usage{
buildScoreCache.bayes(data.df=NULL, data.dists=NULL, group.var=NULL, cor.vars=NULL, 
     dag.banned=NULL, dag.retained=NULL, max.parents=NULL, which.nodes=NULL,
     defn.res=NULL,dry.run=FALSE, max.mode.error=10, verbose=FALSE,
     centre=TRUE, mean=0, prec=0.001, loggam.shape=1,
     loggam.inv.scale=5e-05, max.iters=100, epsabs=1e-7, error.verbose=FALSE,
     epsabs.inner=1e-6,max.iters.inner=100, finite.step.size=1e-7,
     hessian.params=c(1E-04,1E-02), max.iters.hessian=10,
     max.hessian.error=5E-01, factor.brent=1E+02,
     maxiters.hessian.brent=100, num.intervals.brent=100,seed=9062019)
buildScoreCache.mle(data.df=NULL, data.dists=NULL, max.parents=NULL,
     adj.vars=NULL, cor.vars=NULL, dag.banned=NULL, dag.retained=NULL,
     which.nodes=NULL, maxit=100, tol=10^-8, centre=TRUE, defn.res=NULL,
     dry.run=FALSE, verbose=FALSE,seed=9062019)
 
fitAbn.mle(dag=NULL, data.df=NULL, data.dists=NULL, adj.vars=NULL,
      cor.vars=NULL, centre=TRUE, maxit=100, tol=10^-11, verbose=FALSE,seed=9062019)
 
fitAbn.bayes(dag=NULL, data.df=NULL, data.dists=NULL, group.var=NULL,
     cor.vars=NULL,create.graph=FALSE,compute.fixed=FALSE, mean=0, prec=0.001,
     loggam.shape=1,loggam.inv.scale=5e-05,verbose=FALSE, centre=TRUE,
     max.mode.error=10, max.iters=100, epsabs=1e-7, error.verbose=FALSE,
     epsabs.inner=1e-6, max.iters.inner=100, finite.step.size=1E-07,
     hessian.params=c(1E-04,1E-02), max.iters.hessian=10, max.hessian.error=1E-04, 
     factor.brent=1E+02, maxiters.hessian.brent=10, num.intervals.brent=100,
     min.pdf=1E-03,n.grid=100,std.area=TRUE,
     marginal.quantiles=c(0.025,0.25,0.5,0.75,0.975),max.grid.iter=1000,
     marginal.node=NULL, marginal.param=NULL,variate.vec=NULL,seed=9062019)
 
calc.node.inla.glm(child.loc,dag.m.loc,data.df.loc,data.dists.loc,
     ntrials.loc,exposure.loc,compute.fixed.loc,
     mean.intercept.loc,prec.intercept.loc,mean.loc,prec.loc,
     loggam.shape.loc,loggam.inv.scale.loc,verbose.loc)
 
calc.node.inla.glmm(child.loc,dag.m.loc,data.df.loc,data.dists.loc,
     ntrials.loc,exposure.loc,compute.fixed.loc,
     mean.intercept.loc,prec.intercept.loc,mean.loc,prec.loc,
     loggam.shape.loc,loggam.inv.scale.loc,verbose.loc)
          }

\details{
\code{buildScoreCache.mle} and \code{buildScoreCache.bayes} are internal functions called by \code{buildScoreCache}. \cr

\code{fitAbn.mle} and \code{fitAbn.bayes} are internal functions called by \code{fitAbn}. \cr

\code{calc.node.inla.glm} and \code{calc.node.inla.glmm} are internal wrappers to INLA and are called from \code{fitAbn.bayes}.}

\author{Fraser Iain Lewis, Gilles Kratzer}
\keyword{internal}
