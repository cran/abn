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
                      dag.banned=NULL, dag.retained=NULL,  max.parents=NULL,
                      which.nodes=NULL,defn.res=NULL,dry.run=FALSE, centre=TRUE,
                      control = build.control(method = "bayes"), verbose = FALSE)
buildScoreCache.mle(data.df = NULL, data.dists = NULL, max.parents = NULL, 
                    adj.vars = NULL, cor.vars = NULL, dag.banned = NULL, 
                    dag.retained = NULL, which.nodes = NULL, centre = TRUE, 
                    defn.res = NULL, dry.run = FALSE, verbose = FALSE, 
                    control = build.control(method = "mle"))

fitAbn.mle(dag = NULL, data.df = NULL, data.dists = NULL, adj.vars = NULL,
           cor.vars = NULL, centre = TRUE, control = fit.control(method = "mle"),
           verbose = FALSE)
fitAbn.bayes(dag=NULL, data.df=NULL, data.dists=NULL, group.var=NULL,
     cor.vars=NULL, centre=TRUE, compute.fixed=FALSE,
     control = fit.control(method = "bayes"),
     verbose=FALSE)

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
